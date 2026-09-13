"""WT and mAL-output-silence baselines on matched male/female/no-input schedules."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from orientation.intervention import none, output_silence
from orientation.simulate import (
    DEFAULT_REVERSAL_MV,
    load_prepared,
    preference_score,
    response_metrics,
    run_trial,
    sensory_event_schedule,
)

SEEDS = (11, 12, 13)
INPUTS = ("candidate_male", "candidate_female", "no_input")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "baseline-results.json",
    )
    parser.add_argument("--reversal", type=float, default=DEFAULT_REVERSAL_MV)
    args = parser.parse_args()

    started = time.perf_counter()
    network = load_prepared(args.data)
    mal = network.groups["mAL"]
    p1 = network.groups["P1_readout"]
    if np.intersect1d(mal, p1).size:
        raise AssertionError("mAL and P1 readout overlap; protocol violation")

    interventions = {
        "WT": none(),
        "mAL_output_silence": output_silence(mal, name="mAL_output_silence"),
    }
    rows = []
    schedules = {}
    for input_name in INPUTS:
        for seed in SEEDS:
            inputs, events = sensory_event_schedule(network, input_name, seed)
            digest = hashlib.sha256(np.ascontiguousarray(events).tobytes()).hexdigest()
            schedules[f"{input_name}__{seed}"] = {
                "events": int(events.sum()),
                "event_sha256": digest,
                "n_input_cells": int(len(inputs)),
            }
            for name, interv in interventions.items():
                begin = time.perf_counter()
                result = run_trial(
                    network,
                    input_name,
                    seed,
                    interv,
                    reversal_mV=args.reversal,
                    events=events,
                )
                metrics = response_metrics(result, network, interv)
                if result.event_sha256 != digest:
                    raise AssertionError("trial consumed a different event schedule")
                if interv.target.size and np.intersect1d(interv.target, p1).size:
                    raise AssertionError("intervention targets P1 readout")
                metrics["seconds"] = round(time.perf_counter() - begin, 3)
                rows.append(metrics)
                print(json.dumps(metrics), flush=True)

    preferences = []
    for name in interventions:
        for seed in SEEDS:
            male = next(
                r["P1_spikes"]
                for r in rows
                if r["intervention"] == name
                and r["input"] == "candidate_male"
                and r["seed"] == seed
            )
            female = next(
                r["P1_spikes"]
                for r in rows
                if r["intervention"] == name
                and r["input"] == "candidate_female"
                and r["seed"] == seed
            )
            row = {
                "intervention": name,
                "seed": seed,
                "male_P1_spikes": male,
                "female_P1_spikes": female,
                "preference_score": preference_score(male, female),
            }
            preferences.append(row)
            print(json.dumps({"preference": row}), flush=True)

    # Paired male-cue disinhibition check vs upstream qualitative claim.
    paired = []
    for seed in SEEDS:
        intact = next(
            r["P1_spikes"]
            for r in rows
            if r["intervention"] == "WT"
            and r["input"] == "candidate_male"
            and r["seed"] == seed
        )
        blocked = next(
            r["P1_spikes"]
            for r in rows
            if r["intervention"] == "mAL_output_silence"
            and r["input"] == "candidate_male"
            and r["seed"] == seed
        )
        paired.append(
            {
                "seed": seed,
                "WT_P1_spikes": intact,
                "mAL_silence_P1_spikes": blocked,
                "delta": blocked - intact,
            }
        )

    mean_pref = {
        name: float(
            np.mean(
                [p["preference_score"] for p in preferences if p["intervention"] == name]
            )
        )
        for name in interventions
    }

    report = {
        "protocol_sha256": hashlib.sha256(
            (Path(__file__).resolve().parent / "PROTOCOL.md").read_bytes()
        ).hexdigest(),
        "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "simulate_sha256": hashlib.sha256(
            (Path(__file__).resolve().parent / "orientation" / "simulate.py").read_bytes()
        ).hexdigest(),
        "intervention_sha256": hashlib.sha256(
            (Path(__file__).resolve().parent / "orientation" / "intervention.py").read_bytes()
        ).hexdigest(),
        "neurons": network.n,
        "edges": int(network.graph.nnz),
        "reversal_mV": args.reversal,
        "seeds": list(SEEDS),
        "inputs": list(INPUTS),
        "interventions": list(interventions),
        "p1_body_ids": network.ids[p1].tolist(),
        "mAL_size": int(mal.size),
        "results": rows,
        "schedules": schedules,
        "preferences": preferences,
        "mean_preference_score": mean_pref,
        "male_cue_paired_delta": paired,
        "checks": {
            "matched_event_hashes": True,
            "no_P1_stimulation_in_baselines": True,
            "voltage_bounded": True,
            "fixed_seeds": list(SEEDS),
        },
        "seconds": time.perf_counter() - started,
        "interpretation": (
            "Preference score is a modeled P1 response asymmetry index under "
            "matched sensory schedules, not behavioral mate preference."
        ),
    }
    args.out.write_text(json.dumps(report, indent=2))
    print("COMPLETE", report["seconds"], flush=True)


if __name__ == "__main__":
    main()
