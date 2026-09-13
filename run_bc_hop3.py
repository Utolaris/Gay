"""Run frozen S1–S2 block (see BC_HOP3_PROTOCOL.md). Female pathway intact."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np

from orientation.experimental import load_ex_network, run_ex, type_indices
from orientation.simulate import (
    Network,
    apply_sign_overrides,
    load_prepared,
    preference_score,
)

SEEDS = (11, 12, 13)
INPUTS = ("candidate_male", "candidate_female")


def _limit_threads() -> None:
    try:
        from numba import set_num_threads

        set_num_threads(1)
    except Exception:
        pass


_NET_CACHE: dict[str, Network] = {}
_EX_CACHE: dict[str, Network] = {}


def get_net(data_dir: str) -> Network:
    if data_dir not in _NET_CACHE:
        _NET_CACHE[data_dir] = load_prepared(data_dir)
    return _NET_CACHE[data_dir]


def get_ex(data_dir: str) -> Network:
    if data_dir not in _EX_CACHE:
        _EX_CACHE[data_dir] = load_ex_network(Path(data_dir))
    return _EX_CACHE[data_dir]


def _worker(job: dict[str, Any]) -> dict[str, Any]:
    _limit_threads()
    data_dir = job["data_dir"]
    net = get_net(data_dir)
    exnet = get_ex(data_dir)
    signs = apply_sign_overrides(net)
    for tname, sval in job.get("sign_set", []):
        signs[type_indices(net, Path(data_dir), [tname])] = int(sval)

    n = exnet.n
    gain = np.ones(n, dtype=np.float64)
    if job.get("mal_silence"):
        gain[exnet.groups["mAL"]] = 0.0
    for t in job.get("silence_types", []):
        gain[type_indices(exnet, Path(data_dir), [t])] = 0.0
    for t, gval in job.get("gain_types", []):
        gain[type_indices(exnet, Path(data_dir), [t])] = float(gval)
    tonic = np.zeros(n, dtype=np.float64)
    for t, mv in job.get("tonic_types", []):
        tonic[type_indices(exnet, Path(data_dir), [t])] = float(mv)

    metrics = run_ex(
        exnet,
        job["input"],
        job["seed"],
        gain=gain,
        tonic=tonic,
        ei=-70.0,
        use_arousal=False,
        encoding="full",
        signs=signs,
    )
    metrics["block"] = job["block"]
    metrics["condition"] = job["condition"]
    metrics["mal_silence"] = bool(job.get("mal_silence", False))
    metrics["sign_set"] = job.get("sign_set", [])
    return metrics


def build_jobs(data: Path) -> list[dict[str, Any]]:
    jobs: list[dict[str, Any]] = []

    def add(block, condition, **kw):
        for input_name in INPUTS:
            for seed in SEEDS:
                jobs.append(
                    {
                        "data_dir": str(data),
                        "block": block,
                        "condition": condition,
                        "input": input_name,
                        "seed": seed,
                        "sign_set": [],
                        "mal_silence": False,
                        "silence_types": [],
                        "gain_types": [],
                        "tonic_types": [],
                        **kw,
                    }
                )

    variants = {
        "base": [],
        "bc_exc": [("AN09B017b", 1), ("AN09B017c", 1)],
        "bc_zero": [("AN09B017b", 0), ("AN09B017c", 0)],
        "bcg_exc": [
            ("AN09B017b", 1),
            ("AN09B017c", 1),
            ("AN09B017g", 1),
        ],
    }
    for vname, sset in variants.items():
        add("S1", f"{vname}_WT", sign_set=sset)
        add("S1", f"{vname}_mAL", sign_set=sset, mal_silence=True)

    add("S2", "WT")
    add("S2", "mAL", mal_silence=True)
    add("S2", "c_mAL", mal_silence=True, silence_types=["AN09B017c"])
    add(
        "S2",
        "bc_mAL",
        mal_silence=True,
        silence_types=["AN09B017b", "AN09B017c"],
    )
    add(
        "S2",
        "bcd_mAL",
        mal_silence=True,
        silence_types=["AN09B017b", "AN09B017c", "AN09B017d"],
    )
    add("S2", "a008_mAL", mal_silence=True, tonic_types=[("AN03A008", 6.0)])
    add(
        "S2",
        "a008_bc_mAL",
        mal_silence=True,
        silence_types=["AN09B017b", "AN09B017c"],
        tonic_types=[("AN03A008", 6.0)],
    )
    add(
        "S2",
        "a008_c_mAL",
        mal_silence=True,
        silence_types=["AN09B017c"],
        tonic_types=[("AN03A008", 6.0)],
    )
    add("S2", "in11a_mAL", mal_silence=True, silence_types=["IN05B011a"])
    add(
        "S2",
        "in11a_bc_mAL",
        mal_silence=True,
        silence_types=["IN05B011a", "AN09B017b", "AN09B017c"],
    )
    add("S2", "d_g25_mAL", mal_silence=True, gain_types=[("AN09B017d", 0.25)])
    add(
        "S2",
        "a008_in11a_mAL",
        mal_silence=True,
        silence_types=["IN05B011a"],
        tonic_types=[("AN03A008", 6.0)],
    )
    return jobs


def summarize(rows: list[dict]) -> list[dict]:
    by: dict[tuple, list] = {}
    for r in rows:
        by.setdefault((r["block"], r["condition"]), []).append(r)
    summary = []
    for (block, cond), items in by.items():
        def pick(inp, field="P1_spikes"):
            return [
                next(x[field] for x in items if x["input"] == inp and x["seed"] == s)
                for s in SEEDS
            ]

        m, f = pick("candidate_male"), pick("candidate_female")
        scores = [preference_score(a, b) for a, b in zip(m, f)]
        summary.append(
            {
                "block": block,
                "condition": cond,
                "male_P1": m,
                "female_P1": f,
                "male_broad": pick("candidate_male", "broad_P1_spikes"),
                "female_broad": pick("candidate_female", "broad_P1_spikes"),
                "mean_preference": float(np.mean(scores)),
                "male_gt_female_all_seeds": all(a > b for a, b in zip(m, f)),
                "female_bit_identical_to_ref": None,
            }
        )
    refs = {s["condition"]: s for s in summary if s["block"] == "S2"}
    for s in summary:
        if s["block"] != "S2":
            continue
        ref_name = "WT" if s["condition"] == "WT" else "mAL"
        ref = refs.get(ref_name)
        if ref:
            s["female_bit_identical_to_ref"] = s["female_P1"] == ref["female_P1"]
            s["delta_male_vs_ref"] = [a - b for a, b in zip(s["male_P1"], ref["male_P1"])]
            s["delta_female_vs_ref"] = [
                a - b for a, b in zip(s["female_P1"], ref["female_P1"])
            ]
    summary.sort(key=lambda s: (s["block"], s["condition"]))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "bc-hop3-results.json",
    )
    args = parser.parse_args()
    started = time.perf_counter()
    jobs = build_jobs(args.data)
    print(f"jobs={len(jobs)}", flush=True)
    nw = args.workers if args.workers is not None else max(1, min(6, len(jobs)))
    rows: list[dict | None] = [None] * len(jobs)
    with ProcessPoolExecutor(max_workers=nw) as pool:
        futs = {pool.submit(_worker, j): i for i, j in enumerate(jobs)}
        done = 0
        for fut in as_completed(futs):
            i = futs[fut]
            try:
                rows[i] = fut.result()
            except Exception as exc:
                rows[i] = {
                    "block": jobs[i]["block"],
                    "condition": jobs[i]["condition"],
                    "input": jobs[i]["input"],
                    "seed": jobs[i]["seed"],
                    "error": repr(exc),
                }
            done += 1
            if done % 20 == 0:
                print(f"completed {done}/{len(jobs)}", flush=True)
    ok = [r for r in rows if r is not None and "error" not in r]
    summary = summarize(ok) if ok else []
    flipped = [s for s in summary if s["male_gt_female_all_seeds"]]
    report = {
        "kind": "S_block_bc_sign_hop3",
        "protocol_sha256": hashlib.sha256(
            (Path(__file__).resolve().parent / "BC_HOP3_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        "seeds": list(SEEDS),
        "workers": nw,
        "n_jobs": len(jobs),
        "n_ok": len(ok),
        "n_error": len(rows) - len(ok),
        "results": rows,
        "summary": summary,
        "phenotype_C_conditions": [s["condition"] for s in flipped],
        "seconds": time.perf_counter() - started,
    }
    args.out.write_text(json.dumps(report, indent=2))
    print("=== SUMMARY ===", flush=True)
    for s in summary:
        print(
            json.dumps(
                {
                    "block": s["block"],
                    "cond": s["condition"],
                    "M": s["male_P1"],
                    "F": s["female_P1"],
                    "pref": round(s["mean_preference"], 3),
                    "C": s["male_gt_female_all_seeds"],
                    "F_same": s.get("female_bit_identical_to_ref"),
                }
            ),
            flush=True,
        )
    if report["n_error"]:
        print("ERRORS", report["n_error"], flush=True)
        for r in rows:
            if r and "error" in r:
                print(r, flush=True)
    print("COMPLETE", report["seconds"], flush=True)


if __name__ == "__main__":
    main()
