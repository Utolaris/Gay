"""Benchmark trial-level parallel workers: wall time, CPU, peak RSS."""
from __future__ import annotations

import argparse
import json
import os
import resource
import time
from pathlib import Path

from orientation.intervention import none, output_silence
from orientation.parallel import TrialSpec, run_trial_specs
from orientation.simulate import load_prepared

SEEDS = (11, 12, 13)
INPUTS = ("candidate_male", "candidate_female", "no_input")


def peak_rss_mb() -> float:
    # ru_maxrss is KB on Linux, bytes on macOS
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys_platform_is_darwin():
        return rss / (1024 * 1024)
    return rss / 1024


def sys_platform_is_darwin() -> bool:
    import sys

    return sys.platform == "darwin"


def child_peak_rss_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if sys_platform_is_darwin():
        return rss / (1024 * 1024)
    return rss / 1024


def cpu_times() -> dict:
    self_ru = resource.getrusage(resource.RUSAGE_SELF)
    child_ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        "self_user_s": self_ru.ru_utime,
        "self_sys_s": self_ru.ru_stime,
        "children_user_s": child_ru.ru_utime,
        "children_sys_s": child_ru.ru_stime,
        "total_cpu_s": (
            self_ru.ru_utime
            + self_ru.ru_stime
            + child_ru.ru_utime
            + child_ru.ru_stime
        ),
    }


def build_specs(data: Path):
    net = load_prepared(data)
    interventions = {
        "WT": none(),
        "mAL_output_silence": output_silence(net.groups["mAL"], name="mAL_output_silence"),
    }
    specs = []
    for input_name in INPUTS:
        for seed in SEEDS:
            for name, interv in interventions.items():
                specs.append(
                    TrialSpec.from_intervention(data, input_name, seed, interv, label=name)
                )
    return specs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument(
        "--worker-list",
        type=str,
        default="1,2,4,6,8",
        help="comma-separated worker counts",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "parallel-benchmark.json",
    )
    args = parser.parse_args()

    # Warmup: compile Numba + ensure mmap files once
    specs = build_specs(args.data)
    print(f"warmup n_trials={len(specs)}...", flush=True)
    run_trial_specs(specs[:1], workers=1, progress=False)

    results = []
    for w in [int(x) for x in args.worker_list.split(",") if x.strip()]:
        print(f"benchmark workers={w}...", flush=True)
        cpu0 = cpu_times()
        t0 = time.perf_counter()
        run_trial_specs(specs, workers=w, progress=False)
        wall = time.perf_counter() - t0
        cpu1 = cpu_times()
        delta_cpu = cpu1["total_cpu_s"] - cpu0["total_cpu_s"]
        row = {
            "workers": w,
            "n_trials": len(specs),
            "wall_seconds": wall,
            "cpu_seconds": delta_cpu,
            "cpu_utilization": delta_cpu / wall if wall else None,
            "self_peak_rss_mb": peak_rss_mb(),
            "children_peak_rss_mb": child_peak_rss_mb(),
            "cpu_count": os.cpu_count(),
        }
        results.append(row)
        print(json.dumps(row), flush=True)

    # speedup vs workers=1
    base = next((r for r in results if r["workers"] == 1), None)
    for r in results:
        if base and r["wall_seconds"] > 0:
            r["speedup_vs_1"] = base["wall_seconds"] / r["wall_seconds"]

    report = {"results": results, "cpu_count": os.cpu_count()}
    args.out.write_text(json.dumps(report, indent=2))
    print("COMPLETE", flush=True)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
