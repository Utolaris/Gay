"""Compare orientation baselines against committed upstream bounded-route summary."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path)
    parser.add_argument("summary", type=Path)
    args = parser.parse_args()
    base = json.loads(args.baseline.read_text())
    summary = json.loads(args.summary.read_text())
    rows = base["results"]
    pairs = [
        p
        for p in summary["bounded_pairs"]
        if p["input"] in ("candidate_male", "candidate_female") and p["reversal_mV"] == -70
    ]
    mismatches = []
    for p in pairs:
        wt = next(
            r
            for r in rows
            if r["intervention"] == "WT"
            and r["input"] == p["input"]
            and r["seed"] == p["seed"]
        )
        sil = next(
            r
            for r in rows
            if r["intervention"] == "mAL_output_silence"
            and r["input"] == p["input"]
            and r["seed"] == p["seed"]
        )
        if wt["P1_spikes"] != p["intact_spikes"] or sil["P1_spikes"] != p["blocked_spikes"]:
            mismatches.append(
                {
                    "pair": p,
                    "got_intact": wt["P1_spikes"],
                    "got_blocked": sil["P1_spikes"],
                }
            )
    payload = {
        "compared_pairs": len(pairs),
        "mismatches": mismatches,
        "match": not mismatches,
    }
    print(json.dumps(payload, indent=2))
    if mismatches:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
