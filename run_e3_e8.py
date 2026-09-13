"""Run frozen E3–E8 grids (see E3_E8_PROTOCOL.md). Female intact except E8 sensitivity."""
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
from orientation.parallel import default_worker_count
from orientation.simulate import preference_score

SEEDS = (11, 12, 13)
INPUTS = ("candidate_male", "candidate_female")


def _limit_threads() -> None:
    try:
        from numba import set_num_threads

        set_num_threads(1)
    except Exception:
        pass


_NET_CACHE: dict[tuple[str, float], Any] = {}


def _worker(job: dict[str, Any]) -> dict[str, Any]:
    _limit_threads()
    data = Path(job["data_dir"])
    e8 = float(job.get("e8_mult", 1.0))
    cache_key = (str(data.resolve()), e8)
    net = _NET_CACHE.get(cache_key)
    if net is None:
        net = load_ex_network(data, e8_mult=e8)
        _NET_CACHE[cache_key] = net
    mal = net.groups["mAL"]
    n = net.n
    gain = np.ones(n, dtype=np.float64)
    for idxs in job.get("gain0", []):
        gain[np.asarray(idxs, dtype=np.int32)] = 0.0
    for idxs, gval in job.get("gain_set", []):
        gain[np.asarray(idxs, dtype=np.int32)] = gval
    if job.get("mal_silence"):
        gain[mal] = 0.0
    relay = np.zeros(n, dtype=np.bool_)
    for idxs in job.get("relay_idx", []):
        relay[np.asarray(idxs, dtype=np.int32)] = True
    type_maps = job.get("type_maps")
    tmaps = None
    if type_maps:
        tmaps = {k: np.asarray(v, dtype=np.int32) for k, v in type_maps.items()}
    metrics = run_ex(
        net,
        job["input"],
        job["seed"],
        gain=gain,
        ei=-70.0,
        use_arousal=bool(job.get("use_arousal", False)),
        arousal_global=bool(job.get("arousal_global", True)),
        encoding=job.get("encoding", "full"),
        type_maps=tmaps,
    )
    metrics["experiment"] = job["experiment"]
    metrics["condition"] = job["condition"]
    metrics["mal_silence"] = bool(job.get("mal_silence", False))
    metrics["encoding"] = job.get("encoding", "full")
    metrics["e8_mult"] = float(job.get("e8_mult", 1.0))
    return metrics


def build_jobs(data: Path) -> list[dict[str, Any]]:
    # Static type maps from a base network (indices stable across E8 copies)
    net = load_ex_network(data)
    mal = net.groups["mAL"]
    b = type_indices(net, data, ["AN09B017b"])
    c = type_indices(net, data, ["AN09B017c"])
    d = type_indices(net, data, ["AN09B017d"])
    bcd = np.unique(np.concatenate([b, c, d]))
    hub035 = type_indices(net, data, ["AN05B035"])
    hub11a = type_indices(net, data, ["IN05B011a"])
    hub11b = type_indices(net, data, ["IN05B011b"])
    hubs = np.unique(np.concatenate([hub035, hub11a, hub11b]))
    relay = type_indices(
        net, data, ["AN09B017e", "AN09B017f", "AN03A008", "AN09B002"]
    )
    type_maps = {
        "male_a": type_indices(net, data, ["LgLG6"]),
        "female_a": type_indices(net, data, ["LgLG5"]),
        "male_b": type_indices(net, data, ["LgLG7"]),
        "female_b": type_indices(net, data, ["LgLG8"]),
    }
    # JSON-serializable copies for workers
    tmaps_ser = {k: v.tolist() for k, v in type_maps.items()}

    jobs: list[dict[str, Any]] = []

    def add(experiment, condition, **kw):
        for input_name in INPUTS:
            for seed in SEEDS:
                jobs.append(
                    {
                        "data_dir": str(data),
                        "experiment": experiment,
                        "condition": condition,
                        "input": input_name,
                        "seed": seed,
                        "gain0": [],
                        "gain_set": [],
                        "mal_silence": False,
                        "relay_idx": [],
                        "use_arousal": False,
                        "arousal_global": True,
                        "encoding": "full",
                        "e8_mult": 1.0,
                        "type_maps": None,
                        **kw,
                    }
                )

    # --- E3 ---
    add("E3", "WT")
    add("E3", "mAL_ref", mal_silence=True)
    add("E3", "b_sil", gain0=[b.tolist()])
    add("E3", "c_sil", gain0=[c.tolist()])
    add("E3", "d_sil", gain0=[d.tolist()])
    add("E3", "bcd_sil", gain0=[bcd.tolist()])
    add("E3", "bcd_gain0.25", gain_set=[(bcd.tolist(), 0.25)])
    add("E3", "bcd_sil_mAL", gain0=[bcd.tolist()], mal_silence=True)
    add("E3", "bcd_gain0.25_mAL", gain_set=[(bcd.tolist(), 0.25)], mal_silence=True)

    # --- E4 ---
    add("E4", "hub035_sil", gain0=[hub035.tolist()])
    add("E4", "hub11a_sil", gain0=[hub11a.tolist()])
    add("E4", "hub11b_sil", gain0=[hub11b.tolist()])
    add("E4", "hubs_sil", gain0=[hubs.tolist()])
    add("E4", "hubs_gain0.5", gain_set=[(hubs.tolist(), 0.5)])
    add("E4", "hubs_gain0.5_mAL", gain_set=[(hubs.tolist(), 0.5)], mal_silence=True)

    # --- E5 ---
    add("E5", "arousal_global", use_arousal=True, arousal_global=True)
    add(
        "E5",
        "arousal_global_mAL",
        use_arousal=True,
        arousal_global=True,
        mal_silence=True,
    )
    add(
        "E5",
        "arousal_malerelay",
        use_arousal=True,
        arousal_global=False,
        relay_idx=[relay.tolist()],
    )
    add(
        "E5",
        "arousal_malerelay_mAL",
        use_arousal=True,
        arousal_global=False,
        relay_idx=[relay.tolist()],
        mal_silence=True,
    )

    # --- E7 ---
    for enc in ("full", "half_random", "type_only_a", "type_only_b"):
        add("E7", f"{enc}_WT", encoding=enc, type_maps=tmaps_ser)
        add("E7", f"{enc}_mAL", encoding=enc, type_maps=tmaps_ser, mal_silence=True)

    # --- E8 (separate e8_mult; worker loads modified graph) ---
    for mult in (10.0, 100.0):
        add("E8", f"e8_x{int(mult)}_WT", e8_mult=mult)
        add("E8", f"e8_x{int(mult)}_mAL", e8_mult=mult, mal_silence=True)

    return jobs


def summarize(rows: list[dict]) -> list[dict]:
    by: dict[tuple, list] = {}
    for r in rows:
        key = (r["experiment"], r["condition"])
        by.setdefault(key, []).append(r)
    summary = []
    for (ex, cond), items in by.items():
        # order seeds 11,12,13 per input
        def pick(inp):
            return [
                next(
                    x["P1_spikes"]
                    for x in items
                    if x["input"] == inp and x["seed"] == s
                )
                for s in SEEDS
            ]

        m, f = pick("candidate_male"), pick("candidate_female")
        scores = [preference_score(a, b) for a, b in zip(m, f)]
        bm = [
            next(
                x["broad_P1_spikes"]
                for x in items
                if x["input"] == "candidate_male" and x["seed"] == s
            )
            for s in SEEDS
        ]
        bf = [
            next(
                x["broad_P1_spikes"]
                for x in items
                if x["input"] == "candidate_female" and x["seed"] == s
            )
            for s in SEEDS
        ]
        summary.append(
            {
                "experiment": ex,
                "condition": cond,
                "male_P1": m,
                "female_P1": f,
                "male_broad": bm,
                "female_broad": bf,
                "mean_preference": float(np.mean(scores)),
                "male_gt_female_all_seeds": all(a > b for a, b in zip(m, f)),
            }
        )
    summary.sort(key=lambda s: (s["experiment"], s["condition"]))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "e3-e8-results.json",
    )
    args = parser.parse_args()
    started = time.perf_counter()
    jobs = build_jobs(args.data)
    print(f"jobs={len(jobs)}", flush=True)
    nw = default_worker_count(len(jobs)) if args.workers is None else args.workers
    # more trials → allow more workers up to 6
    if args.workers is None:
        nw = max(1, min(6, len(jobs)))
    rows: list[dict | None] = [None] * len(jobs)
    with ProcessPoolExecutor(max_workers=nw) as pool:
        futs = {pool.submit(_worker, j): i for i, j in enumerate(jobs)}
        done = 0
        for fut in as_completed(futs):
            i = futs[fut]
            rows[i] = fut.result()
            done += 1
            if done % 20 == 0:
                print(f"completed {done}/{len(jobs)}", flush=True)
    rows = [r for r in rows if r is not None]
    summary = summarize(rows)
    flipped = [s for s in summary if s["male_gt_female_all_seeds"]]
    report = {
        "kind": "E3_E8_mechanism_grid",
        "protocol_sha256": hashlib.sha256(
            (Path(__file__).resolve().parent / "E3_E8_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        "seeds": list(SEEDS),
        "workers": nw,
        "n_jobs": len(jobs),
        "female_pathway_intact_except_E8": True,
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
                    "ex": s["experiment"],
                    "cond": s["condition"],
                    "M": s["male_P1"],
                    "F": s["female_P1"],
                    "pref": round(s["mean_preference"], 3),
                    "C": s["male_gt_female_all_seeds"],
                }
            ),
            flush=True,
        )
    print("COMPLETE", report["seconds"], flush=True)


if __name__ == "__main__":
    main()
