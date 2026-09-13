"""E2: amplify male-specific excitatory relays with female pathway intact.

See E2_PROTOCOL.md. Uses per-neuron output gain + tonic vector so mAL silence
and relay activation can be combined without changing baseline simulate.py.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyarrow.feather as feather
from numba import njit
from scipy.sparse import load_npz

from orientation.simulate import (
    BIN_MS,
    DELAY_MS,
    DT_MS,
    EXTERNAL_MV,
    MEMBRANE_TAU_MS,
    REFRACTORY_MS,
    REST_MV,
    STIM_END_MS,
    STIM_START_MS,
    SYN_TAU_MS,
    THRESH_MV,
    WEIGHT_SCALE,
    DURATION_MS,
    load_prepared,
    preference_score,
    sensory_event_schedule,
)

SEEDS = (11, 12, 13)
INPUTS = ("candidate_male", "candidate_female")
REVERSAL = -70.0
TONIC_MV = 6.0

RELAY_SETS = {
    "vAB3ef": ["AN09B017e", "AN09B017f"],
    "AN03A008": ["AN03A008"],
    "AN09B002": ["AN09B002"],
    "all4": ["AN09B017e", "AN09B017f", "AN03A008", "AN09B002"],
}


@njit(cache=True)
def _sim(
    indptr,
    dest,
    weights,
    signs,
    sensory_inputs,
    sensory_events,
    neuron_gain,
    tonic,
    ei,
):
    n = signs.shape[0]
    steps = sensory_events.shape[0]
    delay = int(round(DELAY_MS / DT_MS))
    ref = int(round(REFRACTORY_MS / DT_MS))
    v = np.full(n, REST_MV)
    ge = np.zeros(n)
    hi = np.zeros(n)
    release = np.zeros(n, np.int32)
    pending_e = np.zeros((delay + 1, n))
    pending_i = np.zeros((delay + 1, n))
    stim_mask = np.zeros(n, np.bool_)
    for j in range(sensory_inputs.shape[0]):
        stim_mask[sensory_inputs[j]] = True
    n_bins = int(np.ceil(steps * DT_MS / BIN_MS))
    counts = np.zeros((n_bins, n), np.int32)
    eg = np.exp(-DT_MS / SYN_TAU_MS)
    for step in range(steps):
        t = step * DT_MS
        slot = step % (delay + 1)
        for i in range(n):
            inc_e = pending_e[slot, i]
            inc_i = pending_i[slot, i]
            pending_e[slot, i] = 0.0
            pending_i[slot, i] = 0.0
            if step >= release[i]:
                ge[i] += inc_e
                hi[i] += inc_i
                rate = (1.0 + hi[i]) / MEMBRANE_TAU_MS
                decay = np.exp(-rate * DT_MS)
                equilibrium = (REST_MV + hi[i] * ei + tonic[i]) / (1.0 + hi[i])
                denom = rate - 1.0 / SYN_TAU_MS
                if abs(denom) > 1e-10:
                    coupling = (eg - decay) / (MEMBRANE_TAU_MS * denom)
                else:
                    coupling = DT_MS * decay / MEMBRANE_TAU_MS
                v[i] = equilibrium + (v[i] - equilibrium) * decay + ge[i] * coupling
                ge[i] *= eg
                hi[i] *= eg
        for j in range(sensory_inputs.shape[0]):
            if sensory_events[step, j]:
                v[sensory_inputs[j]] += EXTERNAL_MV
        for i in range(n):
            if step >= release[i] and v[i] > THRESH_MV:
                bin_i = int(t / BIN_MS)
                if bin_i >= n_bins:
                    bin_i = n_bins - 1
                counts[bin_i, i] += 1
                v[i] = REST_MV
                ge[i] = 0.0
                hi[i] = 0.0
                release[i] = step + (0 if stim_mask[i] else ref)
                if signs[i] == 0 or neuron_gain[i] == 0.0:
                    continue
                slot2 = (step + delay) % (delay + 1)
                g = neuron_gain[i]
                for e in range(indptr[i], indptr[i + 1]):
                    target = dest[e]
                    w = weights[e] * WEIGHT_SCALE * g
                    if signs[i] > 0:
                        pending_e[slot2, target] += w
                    else:
                        pending_i[slot2, target] += w / (REST_MV - ei)
    return counts


def type_indices(net, data_dir: Path, types: list[str]) -> np.ndarray:
    ann = feather.read_table(data_dir / "annotations.feather").to_pylist()
    by_id = {r["bodyId"]: r for r in ann}
    want = set(types)
    idxs = [
        i
        for i, bid in enumerate(net.ids)
        if (by_id.get(int(bid)) or {}).get("type") in want
    ]
    return np.asarray(idxs, dtype=np.int32)


@dataclass(frozen=True)
class Cond:
    name: str
    mal_silence: bool
    relay_types: tuple[str, ...] = ()
    tonic_mV: float = 0.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "e2-results.json",
    )
    args = parser.parse_args()
    data_dir = args.data
    started = time.perf_counter()

    net = load_prepared(data_dir)
    # baseline signs (female + vAB3 excitatory)
    signs = net.signs.copy()
    signs[net.groups["vAB3"]] = 1
    signs[net.groups["candidate_female"]] = 1
    p1 = net.groups["P1_readout"]
    mal = net.groups["mAL"]
    if np.intersect1d(mal, p1).size or np.intersect1d(
        type_indices(net, data_dir, ["AN09B017e", "AN09B017f"]), p1
    ).size:
        raise AssertionError("intervention target overlaps P1")

    conds = [
        Cond("WT", mal_silence=False),
        Cond("mAL_silence", mal_silence=True),
    ]
    for name, types in RELAY_SETS.items():
        conds.append(Cond(f"{name}_tonic6", mal_silence=False, relay_types=tuple(types), tonic_mV=TONIC_MV))
        conds.append(
            Cond(
                f"{name}_tonic6_mAL",
                mal_silence=True,
                relay_types=tuple(types),
                tonic_mV=TONIC_MV,
            )
        )

    relay_cache = {name: type_indices(net, data_dir, types) for name, types in RELAY_SETS.items()}
    for name, idxs in relay_cache.items():
        if idxs.size == 0:
            raise AssertionError(f"empty relay set {name}")
        print(f"relay {name}: n={idxs.size}", flush=True)

    rows = []
    for cond in conds:
        gain = np.ones(net.n, dtype=np.float64)
        if cond.mal_silence:
            gain[mal] = 0.0
        tonic = np.zeros(net.n, dtype=np.float64)
        if cond.relay_types:
            # map condensed name back
            key = {
                "vAB3ef": "vAB3ef",
                "AN03A008": "AN03A008",
                "AN09B002": "AN09B002",
                "all4": "all4",
            }[cond.name.split("_tonic")[0]]
            tonic[relay_cache[key]] = cond.tonic_mV
        for input_name in INPUTS:
            for seed in SEEDS:
                begin = time.perf_counter()
                inputs, events = sensory_event_schedule(net, input_name, seed)
                digest = hashlib.sha256(np.ascontiguousarray(events).tobytes()).hexdigest()
                counts = _sim(
                    net.graph.indptr,
                    net.graph.indices,
                    net.graph.data,
                    signs,
                    inputs,
                    events,
                    gain,
                    tonic,
                    REVERSAL,
                )
                window = counts[5:]
                p1_spikes = int(window[:, p1].sum())
                row = {
                    "condition": cond.name,
                    "mal_silence": cond.mal_silence,
                    "relay_types": list(cond.relay_types),
                    "tonic_mV": cond.tonic_mV,
                    "input": input_name,
                    "seed": seed,
                    "events": int(events.sum()),
                    "event_sha256": digest,
                    "P1_spikes": p1_spikes,
                    "broad_P1_spikes": int(window[:, net.groups["P1_related"]].sum()),
                    "mAL_spikes": int(window[:, mal].sum()),
                    "network_spikes": int(counts.sum()),
                    "active_neurons": int(np.any(counts, axis=0).sum()),
                    "seconds": round(time.perf_counter() - begin, 3),
                }
                rows.append(row)
                print(json.dumps(row), flush=True)

    def p1_of(cond: str, input_name: str, seed: int) -> int:
        return next(
            r["P1_spikes"]
            for r in rows
            if r["condition"] == cond and r["input"] == input_name and r["seed"] == seed
        )

    summary = []
    for cond in conds:
        if cond.name == "WT":
            continue
        male = [p1_of(cond.name, "candidate_male", s) for s in SEEDS]
        female = [p1_of(cond.name, "candidate_female", s) for s in SEEDS]
        wt_m = [p1_of("WT", "candidate_male", s) for s in SEEDS]
        wt_f = [p1_of("WT", "candidate_female", s) for s in SEEDS]
        mal_m = [p1_of("mAL_silence", "candidate_male", s) for s in SEEDS]
        mal_f = [p1_of("mAL_silence", "candidate_female", s) for s in SEEDS]
        scores = [preference_score(m, f) for m, f in zip(male, female)]
        summary.append(
            {
                "condition": cond.name,
                "male_P1": male,
                "female_P1": female,
                "delta_male_vs_WT": [m - w for m, w in zip(male, wt_m)],
                "delta_female_vs_WT": [f - w for f, w in zip(female, wt_f)],
                "delta_male_vs_mAL": [m - a for m, a in zip(male, mal_m)],
                "delta_female_vs_mAL": [f - a for f, a in zip(female, mal_f)],
                "mean_preference": float(np.mean(scores)),
                "male_gt_female_all_seeds": all(m > f for m, f in zip(male, female)),
                "female_unchanged_vs_mAL_all_seeds": all(
                    f == a for f, a in zip(female, mal_f)
                )
                if cond.mal_silence
                else all(f == w for f, w in zip(female, wt_f)),
            }
        )

    report = {
        "kind": "E2_male_excitatory_relay_amplification",
        "protocol_sha256": hashlib.sha256(
            (Path(__file__).resolve().parent / "E2_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "tonic_mV": TONIC_MV,
        "reversal_mV": REVERSAL,
        "seeds": list(SEEDS),
        "female_pathway_intact": True,
        "relay_sizes": {k: int(v.size) for k, v in relay_cache.items()},
        "results": rows,
        "summary": summary,
        "phenotype_C_conditions": [
            s["condition"] for s in summary if s["male_gt_female_all_seeds"]
        ],
        "seconds": time.perf_counter() - started,
    }
    args.out.write_text(json.dumps(report, indent=2))
    print("=== E2 SUMMARY ===", flush=True)
    for s in summary:
        print(
            json.dumps(
                {
                    "cond": s["condition"],
                    "M": s["male_P1"],
                    "F": s["female_P1"],
                    "pref": s["mean_preference"],
                    "C": s["male_gt_female_all_seeds"],
                    "F_stable": s["female_unchanged_vs_mAL_all_seeds"],
                }
            ),
            flush=True,
        )
    print("COMPLETE", report["seconds"], flush=True)


if __name__ == "__main__":
    main()
