"""E1: single-class P1 brake lesions under matched male/female drive.

See E1_PROTOCOL.md. Female pathway and sign map stay intact.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pyarrow.feather as feather

from orientation.intervention import none, output_silence
from orientation.simulate import (
    DEFAULT_REVERSAL_MV,
    load_prepared,
    preference_score,
    run_trial,
    sensory_event_schedule,
)

SEEDS = (11, 12, 13)
INPUTS = ("candidate_male", "candidate_female")

# Frozen in E1_PROTOCOL.md before outcomes.
LESION_SPECS = [
    ("WT", None),
    ("mAL_all", "group:mAL"),
    ("mAL_m8", "type:mAL_m8"),
    ("mAL_m5b", "type:mAL_m5b"),
    ("mAL_m1", "type:mAL_m1"),
    ("mAL_m2b", "type:mAL_m2b"),
    ("mAL_m5c", "type:mAL_m5c"),
    ("SIP100m", "type:SIP100m"),
    ("VES022", "type:VES022"),
    ("SIP112m", "type:SIP112m"),
    ("LH004m", "type:LH004m"),
    ("SIP113m", "type:SIP113m"),
    ("AN09B017b", "type:AN09B017b"),
    ("AN09B017c", "type:AN09B017c"),
    ("AN09B017d", "type:AN09B017d"),
]


def resolve_indices(net, data_dir: Path, spec: str | None) -> np.ndarray:
    if spec is None:
        return np.zeros(0, dtype=np.int32)
    kind, name = spec.split(":", 1)
    if kind == "group":
        return np.asarray(net.groups[name], dtype=np.int32)
    if kind == "type":
        ann = feather.read_table(data_dir / "annotations.feather").to_pylist()
        by_id = {r["bodyId"]: r for r in ann}
        idxs = [
            i
            for i, bid in enumerate(net.ids)
            if (by_id.get(int(bid)) or {}).get("type") == name
        ]
        return np.asarray(idxs, dtype=np.int32)
    raise ValueError(spec)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "e1-results.json",
    )
    parser.add_argument("--reversal", type=float, default=DEFAULT_REVERSAL_MV)
    args = parser.parse_args()

    started = time.perf_counter()
    net = load_prepared(args.data)
    p1 = net.groups["P1_readout"]
    if np.intersect1d(net.groups["mAL"], p1).size:
        raise AssertionError("mAL overlaps P1 readout")

    lesions = []
    for name, spec in LESION_SPECS:
        idxs = resolve_indices(net, args.data, spec)
        if name != "WT" and idxs.size == 0:
            raise AssertionError(f"empty lesion class: {name}")
        if np.intersect1d(idxs, p1).size:
            raise AssertionError(f"lesion {name} targets P1 readout")
        interv = none() if name == "WT" else output_silence(idxs, name=name)
        lesions.append((name, interv, idxs))

    rows = []
    # Pre-generate matched schedules once per (input, seed).
    schedules = {}
    for input_name in INPUTS:
        for seed in SEEDS:
            _, events = sensory_event_schedule(net, input_name, seed)
            digest = hashlib.sha256(np.ascontiguousarray(events).tobytes()).hexdigest()
            schedules[f"{input_name}__{seed}"] = {"events": int(events.sum()), "sha256": digest}

    for name, interv, idxs in lesions:
        for input_name in INPUTS:
            for seed in SEEDS:
                begin = time.perf_counter()
                key = f"{input_name}__{seed}"
                result = run_trial(
                    net,
                    input_name,
                    seed,
                    interv,
                    reversal_mV=args.reversal,
                    events=sensory_event_schedule(net, input_name, seed)[1],
                )
                window = result.counts[5:]
                p1_spikes = int(window[:, p1].sum())
                row = {
                    "lesion": name,
                    "n_cells": int(idxs.size),
                    "input": input_name,
                    "seed": seed,
                    "events": result.events,
                    "event_sha256": result.event_sha256,
                    "P1_spikes": p1_spikes,
                    "broad_P1_spikes": int(window[:, net.groups["P1_related"]].sum()),
                    "mAL_spikes": int(window[:, net.groups["mAL"]].sum()),
                    "network_spikes": int(result.counts.sum()),
                    "active_neurons": int(np.any(result.counts, axis=0).sum()),
                    "voltage_min_mV": result.voltage_min_mV,
                    "seconds": round(time.perf_counter() - begin, 3),
                }
                if result.event_sha256 != schedules[key]["sha256"]:
                    raise AssertionError("schedule mismatch")
                rows.append(row)
                print(json.dumps(row), flush=True)

    def p1_of(lesion: str, input_name: str, seed: int) -> int:
        return next(
            r["P1_spikes"]
            for r in rows
            if r["lesion"] == lesion and r["input"] == input_name and r["seed"] == seed
        )

    contrasts = []
    for name, _, idxs in lesions:
        if name == "WT":
            continue
        for seed in SEEDS:
            dm = p1_of(name, "candidate_male", seed) - p1_of("WT", "candidate_male", seed)
            df = p1_of(name, "candidate_female", seed) - p1_of("WT", "candidate_female", seed)
            m = p1_of(name, "candidate_male", seed)
            f = p1_of(name, "candidate_female", seed)
            contrasts.append(
                {
                    "lesion": name,
                    "n_cells": int(idxs.size),
                    "seed": seed,
                    "male_P1": m,
                    "female_P1": f,
                    "delta_male": dm,
                    "delta_female": df,
                    "delta_male_minus_female": dm - df,
                    "preference_score": preference_score(m, f),
                }
            )

    summary = []
    names = [n for n, _, _ in lesions if n != "WT"]
    for name in names:
        cs = [c for c in contrasts if c["lesion"] == name]
        dm = [c["delta_male"] for c in cs]
        df = [c["delta_female"] for c in cs]
        summary.append(
            {
                "lesion": name,
                "n_cells": cs[0]["n_cells"],
                "male_P1": [c["male_P1"] for c in cs],
                "female_P1": [c["female_P1"] for c in cs],
                "delta_male": dm,
                "delta_female": df,
                "mean_delta_male": float(np.mean(dm)),
                "mean_delta_female": float(np.mean(df)),
                "mean_delta_male_minus_female": float(np.mean(np.array(dm) - np.array(df))),
                "raises_male_all_seeds": all(d > 0 for d in dm),
                "male_specific_all_seeds": all(d > 0 and d > df[i] for i, d in enumerate(dm)),
            }
        )
    summary.sort(key=lambda s: -s["mean_delta_male_minus_female"])

    report = {
        "kind": "E1_P1_brake_lesion_map",
        "protocol_sha256": hashlib.sha256(
            (Path(__file__).resolve().parent / "E1_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reversal_mV": args.reversal,
        "seeds": list(SEEDS),
        "inputs": list(INPUTS),
        "female_pathway_intact": True,
        "neurons": net.n,
        "results": rows,
        "contrasts": contrasts,
        "summary": summary,
        "seconds": time.perf_counter() - started,
    }
    args.out.write_text(json.dumps(report, indent=2))
    print("=== E1 SUMMARY (Δmale−Δfemale desc) ===", flush=True)
    for s in summary:
        print(
            json.dumps(
                {
                    "lesion": s["lesion"],
                    "n": s["n_cells"],
                    "dM": s["mean_delta_male"],
                    "dF": s["mean_delta_female"],
                    "dM-dF": s["mean_delta_male_minus_female"],
                    "male_up_all": s["raises_male_all_seeds"],
                    "male_spec": s["male_specific_all_seeds"],
                }
            ),
            flush=True,
        )
    print("COMPLETE", report["seconds"], flush=True)


if __name__ == "__main__":
    main()
