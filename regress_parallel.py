"""Strict serial vs parallel regression on all pre-registered baseline trials.

Compares every metric except wall-clock ``seconds``. Exit code 1 on mismatch.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from orientation.intervention import none, output_silence
from orientation.parallel import TrialSpec, run_trial_specs
from orientation.simulate import load_prepared

SEEDS = (11, 12, 13)
INPUTS = ("candidate_male", "candidate_female", "no_input")
IGNORE_KEYS = {"seconds"}


def metrics_key(row: dict) -> tuple:
    return (row["intervention"], row["input"], row["seed"])


def comparable(row: dict) -> dict:
    return {k: v for k, v in row.items() if k not in IGNORE_KEYS and k != "label"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "parallel-regression.json",
    )
    args = parser.parse_args()

    net = load_prepared(args.data)
    interventions = {
        "WT": none(),
        "mAL_output_silence": output_silence(net.groups["mAL"], name="mAL_output_silence"),
    }
    specs = []
    for input_name in INPUTS:
        for seed in SEEDS:
            for name, interv in interventions.items():
                specs.append(
                    TrialSpec.from_intervention(
                        args.data, input_name, seed, interv, label=name
                    )
                )

    print(f"trials={len(specs)} serial...", flush=True)
    t0 = time.perf_counter()
    serial = run_trial_specs(specs, workers=1, progress=False)
    t_serial = time.perf_counter() - t0

    print(f"parallel workers={args.workers}...", flush=True)
    t0 = time.perf_counter()
    parallel = run_trial_specs(specs, workers=args.workers, progress=False)
    t_parallel = time.perf_counter() - t0

    # normalize intervention field from label if needed
    for rows in (serial, parallel):
        for r in rows:
            r["intervention"] = r.get("label") or r["intervention"]

    s_map = {metrics_key(r): comparable(r) for r in serial}
    p_map = {metrics_key(r): comparable(r) for r in parallel}
    mismatches = []
    for k in sorted(s_map.keys() | p_map.keys()):
        if k not in s_map or k not in p_map:
            mismatches.append({"key": k, "error": "missing side"})
            continue
        if s_map[k] != p_map[k]:
            diff = {
                field: {"serial": s_map[k].get(field), "parallel": p_map[k].get(field)}
                for field in sorted(set(s_map[k]) | set(p_map[k]))
                if s_map[k].get(field) != p_map[k].get(field)
            }
            mismatches.append({"key": k, "diff": diff})

    report = {
        "n_trials": len(specs),
        "workers": args.workers,
        "serial_seconds": t_serial,
        "parallel_seconds": t_parallel,
        "speedup": t_serial / t_parallel if t_parallel else None,
        "match": not mismatches,
        "mismatches": mismatches,
    }
    args.out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)
    if mismatches:
        print(f"FAIL: {len(mismatches)} mismatches", flush=True)
        sys.exit(1)
    print("PASS: serial and parallel baselines identical", flush=True)


if __name__ == "__main__":
    main()
