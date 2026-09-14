"""OA / VES022 / SIP106m gate on male→P1 (see OA_GATE_PROTOCOL.md)."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from orientation.experimental import load_ex_network, run_ex, type_indices
from orientation.simulate import preference_score

SEEDS = (11, 12, 13)
INPUTS = ("candidate_male", "candidate_female")
HZ = 25.0
CONDS = ("WT", "OA", "ves022", "sip106")


def _limit_threads() -> None:
    try:
        from numba import set_num_threads

        set_num_threads(1)
    except Exception:
        pass


def _worker(job: dict) -> dict:
    _limit_threads()
    data = Path(job["data_dir"])
    net = load_ex_network(data)
    oa = type_indices(net, data, ["OA-", "OAN"])  # prefix won't match; fix below
    # type_indices uses exact type match — expand OA types explicitly
    oa_types = [
        "OA-AL2i1",
        "OA-AL2i2",
        "OA-AL2i3",
        "OA-AL2i4",
        "OA-ASM1",
        "OA-ASM2",
        "OA-ASM3",
        "OA-VPM3",
        "OA-VPM4",
        "OA-VUMa1",
        "OA-VUMa2",
        "OA-VUMa3",
        "OA-VUMa4",
        "OA-VUMa5",
        "OA-VUMa6",
        "OA-VUMa8",
    ]
    oa = type_indices(net, data, oa_types)
    ves = type_indices(net, data, ["VES022"])
    sip = type_indices(net, data, ["SIP106m"])

    n = net.n
    tonic = np.zeros(n, dtype=np.float64)
    # run_ex activation is via tonic array + we need poisson; experimental
    # simulate_ex has no poisson activation. Use tonic 12 mV on targets
    # (crosses rest–threshold gap of 7 mV) — pre-registered as the drive
    # level in this protocol file before outcomes: 12 mV tonic.
    DRIVE_MV = 12.0
    if job["condition"] == "OA":
        tonic[oa] = DRIVE_MV
    elif job["condition"] == "ves022":
        tonic[ves] = DRIVE_MV
    elif job["condition"] == "sip106":
        tonic[sip] = DRIVE_MV

    m = run_ex(
        net,
        job["input"],
        job["seed"],
        gain=None,
        tonic=tonic,
        ei=-70.0,
        use_arousal=False,
    )
    m["condition"] = job["condition"]
    m["drive_mV"] = DRIVE_MV
    m["n_OA"] = int(oa.size)
    m["n_VES022"] = int(ves.size)
    m["n_SIP106m"] = int(sip.size)
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data", type=Path)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out", type=Path, default=Path("oa-gate-results.json"))
    args = ap.parse_args()
    started = time.perf_counter()
    jobs = [
        {
            "data_dir": str(args.data),
            "condition": c,
            "input": inp,
            "seed": s,
        }
        for c in CONDS
        for inp in INPUTS
        for s in SEEDS
    ]
    print(f"jobs={len(jobs)} drive=12mV tonic", flush=True)
    rows: list[dict | None] = [None] * len(jobs)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(_worker, j): i for i, j in enumerate(jobs)}
        done = 0
        for fut in as_completed(futs):
            i = futs[fut]
            try:
                rows[i] = fut.result()
            except Exception as exc:
                rows[i] = {
                    "condition": jobs[i]["condition"],
                    "input": jobs[i]["input"],
                    "seed": jobs[i]["seed"],
                    "error": repr(exc),
                }
            done += 1
            if done % 6 == 0:
                print(f"completed {done}/{len(jobs)}", flush=True)
    ok = [r for r in rows if r is not None and "error" not in r]
    by = {}
    for r in ok:
        by.setdefault(r["condition"], []).append(r)
    summary = []
    for c in CONDS:
        items = sorted(by.get(c, []), key=lambda r: (r["input"], r["seed"]))
        def pick(inp):
            return [
                next(x["P1_spikes"] for x in items if x["input"] == inp and x["seed"] == s)
                for s in SEEDS
            ]
        m, f = pick("candidate_male"), pick("candidate_female")
        scores = [preference_score(a, b) for a, b in zip(m, f)]
        summary.append(
            {
                "condition": c,
                "male_P1": m,
                "female_P1": f,
                "mean_pref": float(np.mean(scores)),
                "C": all(a > b for a, b in zip(m, f)),
            }
        )
    report = {
        "kind": "OA_gate",
        "protocol_sha256": hashlib.sha256(Path("OA_GATE_PROTOCOL.md").read_bytes()).hexdigest(),
        "drive_mV": 12.0,
        "seeds": list(SEEDS),
        "n_ok": len(ok),
        "n_error": len(rows) - len(ok),
        "results": rows,
        "summary": summary,
        "seconds": time.perf_counter() - started,
    }
    args.out.write_text(json.dumps(report, indent=2))
    for s in summary:
        print(json.dumps(s), flush=True)
    print("COMPLETE", report["seconds"], flush=True)


if __name__ == "__main__":
    main()
