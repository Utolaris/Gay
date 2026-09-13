"""POST-HOC exploratory parameter search: conditions that make male P1 > female P1.

This is NOT the pre-registered baseline. It does not replace PROTOCOL.md or
baseline-results.json. Every condition is reported, including failures.

Axes explored (all with matched seeds 11/12/13, reversal -70 mV):
  - sensory event rate asymmetry
  - remove candidate_female excitatory sign override
  - output_gain on candidate_female
  - tonic activation on candidate_male (not P1)
  - combinations with mAL output silence
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numba import njit
from scipy.sparse import load_npz

from orientation.simulate import (
    BIN_MS,
    DELAY_MS,
    DT_MS,
    EXTERNAL_MV,
    MEMBRANE_TAU_MS,
    P1_BODY_IDS,
    REFRACTORY_MS,
    REST_MV,
    STIM_END_MS,
    STIM_START_MS,
    SYN_TAU_MS,
    THRESH_MV,
    VAB3_BODY_IDS,
    WEIGHT_SCALE,
    DURATION_MS,
    preference_score,
)
from orientation.intervention import none, output_gain, output_silence, activate_tonic

SEEDS = (11, 12, 13)
REVERSAL = -70.0


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


def events_for(net, input_name: str, seed: int, total_hz: float) -> tuple[np.ndarray, np.ndarray]:
    steps = int(round(DURATION_MS / DT_MS))
    s0 = int(round(STIM_START_MS / DT_MS))
    s1 = int(round(STIM_END_MS / DT_MS))
    if input_name == "no_input":
        return np.zeros(0, np.int32), np.zeros((steps, 0), dtype=np.bool_)
    inputs = net.groups[input_name]
    rate = total_hz / len(inputs)
    ev = np.zeros((steps, len(inputs)), dtype=np.bool_)
    ev[s0:s1] = np.random.RandomState(seed).random_sample((s1 - s0, len(inputs))) < (
        rate * DT_MS / 1000.0
    )
    return inputs, ev


def base_signs(net, female_excitatory: bool = True) -> np.ndarray:
    signs = net.signs.copy()
    signs[net.groups["vAB3"]] = 1
    if female_excitatory:
        signs[net.groups["candidate_female"]] = 1
    return signs


def gain_vector(net, pairs: dict[str, float]) -> np.ndarray:
    g = np.ones(net.n, dtype=np.float64)
    for name, value in pairs.items():
        g[net.groups[name]] = value
    return g


def tonic_vector(net, pairs: dict[str, float]) -> np.ndarray:
    t = np.zeros(net.n, dtype=np.float64)
    for name, value in pairs.items():
        t[net.groups[name]] = value
    return t


@dataclass
class Cond:
    name: str
    male_hz: float = 5550.0
    female_hz: float = 5550.0
    female_excitatory: bool = True
    gains: dict | None = None
    tonics: dict | None = None

    def describe(self) -> dict:
        return {
            "male_hz": self.male_hz,
            "female_hz": self.female_hz,
            "female_excitatory": self.female_excitatory,
            "gains": self.gains or {},
            "tonics_mV": self.tonics or {},
        }


def run_cond(net, cond: Cond, input_name: str, seed: int) -> dict:
    total = cond.male_hz if input_name == "candidate_male" else cond.female_hz
    inputs, events = events_for(net, input_name, seed, total)
    signs = base_signs(net, cond.female_excitatory)
    gains = gain_vector(net, cond.gains or {})
    tonic = tonic_vector(net, cond.tonics or {})
    counts = _sim(
        net.graph.indptr,
        net.graph.indices,
        net.graph.data,
        signs,
        inputs,
        events,
        gains,
        tonic,
        REVERSAL,
    )
    p1 = net.groups["P1_readout"]
    window = counts[5:]
    p1_spikes = int(window[:, p1].sum())
    return {
        "condition": cond.name,
        "input": input_name,
        "seed": seed,
        "events": int(events.sum()),
        "P1_spikes": p1_spikes,
        "network_spikes": int(counts.sum()),
        "active_neurons": int(np.any(counts, axis=0).sum()),
    }


def main() -> None:
    from orientation.simulate import load_prepared

    root = Path(".experiment-data")
    net = load_prepared(root)
    mal = "mAL"
    female = "candidate_female"
    male = "candidate_male"

    conds = [
        Cond("ref_WT"),
        Cond("ref_mAL_silence", gains={mal: 0.0}),
        # remove female excitatory override
        Cond("no_fem_exc_WT", female_excitatory=False),
        Cond("no_fem_exc_mAL", female_excitatory=False, gains={mal: 0.0}),
        # male sensory rate up
        Cond("male_x2_WT", male_hz=11100.0),
        Cond("male_x2_mAL", male_hz=11100.0, gains={mal: 0.0}),
        Cond("male_x4_WT", male_hz=22200.0),
        Cond("male_x4_mAL", male_hz=22200.0, gains={mal: 0.0}),
        Cond("male_x8_WT", male_hz=44400.0),
        Cond("male_x8_mAL", male_hz=44400.0, gains={mal: 0.0}),
        # female sensory rate down
        Cond("fem_half_WT", female_hz=2775.0),
        Cond("fem_half_mAL", female_hz=2775.0, gains={mal: 0.0}),
        Cond("fem_quarter_WT", female_hz=1387.5),
        Cond("fem_quarter_mAL", female_hz=1387.5, gains={mal: 0.0}),
        # female output gain
        Cond("fem_gain0.5_WT", gains={female: 0.5}),
        Cond("fem_gain0.5_mAL", gains={mal: 0.0, female: 0.5}),
        Cond("fem_gain0.25_mAL", gains={mal: 0.0, female: 0.25}),
        Cond("fem_gain0.0_mAL", gains={mal: 0.0, female: 0.0}),
        # male class tonic activation (not P1)
        Cond("male_tonic5_WT", tonics={male: 5.0}),
        Cond("male_tonic10_WT", tonics={male: 10.0}),
        Cond("male_tonic10_mAL", tonics={male: 10.0}, gains={mal: 0.0}),
        Cond("male_tonic14_mAL", tonics={male: 14.0}, gains={mal: 0.0}),
        # combined: male up + female down + mAL silence
        Cond(
            "combo_male_x2_fem_half_mAL",
            male_hz=11100.0,
            female_hz=2775.0,
            gains={mal: 0.0},
        ),
        Cond(
            "combo_male_x2_fem_gain0.25_mAL",
            male_hz=11100.0,
            gains={mal: 0.0, female: 0.25},
        ),
        Cond(
            "combo_noexc_male_x2_mAL",
            male_hz=11100.0,
            female_excitatory=False,
            gains={mal: 0.0},
        ),
        Cond(
            "combo_tonic10_fem_gain0.25_mAL",
            tonics={male: 10.0},
            gains={mal: 0.0, female: 0.25},
        ),
    ]

    started = time.perf_counter()
    rows = []
    for cond in conds:
        for seed in SEEDS:
            male_r = run_cond(net, cond, "candidate_male", seed)
            female_r = run_cond(net, cond, "candidate_female", seed)
            score = preference_score(male_r["P1_spikes"], female_r["P1_spikes"])
            row = {
                **cond.describe(),
                "condition": cond.name,
                "seed": seed,
                "male_P1_spikes": male_r["P1_spikes"],
                "female_P1_spikes": female_r["P1_spikes"],
                "male_events": male_r["events"],
                "female_events": female_r["events"],
                "male_network_spikes": male_r["network_spikes"],
                "female_network_spikes": female_r["network_spikes"],
                "preference_score": score,
                "male_gt_female": male_r["P1_spikes"] > female_r["P1_spikes"],
            }
            rows.append(row)
            print(json.dumps(row), flush=True)

    # summarize by condition
    summary = []
    by = {}
    for r in rows:
        by.setdefault(r["condition"], []).append(r)
    for name, items in by.items():
        m = [x["male_P1_spikes"] for x in items]
        f = [x["female_P1_spikes"] for x in items]
        s = [x["preference_score"] for x in items]
        wins = sum(1 for x in items if x["male_gt_female"])
        summary.append(
            {
                "condition": name,
                "params": {k: items[0][k] for k in ["male_hz", "female_hz", "female_excitatory", "gains", "tonics_mV"]},
                "male_P1_spikes": m,
                "female_P1_spikes": f,
                "mean_male": float(np.mean(m)),
                "mean_female": float(np.mean(f)),
                "mean_preference": float(np.mean(s)),
                "seeds_male_gt_female": wins,
                "all_seeds_male_gt_female": wins == len(SEEDS),
            }
        )
    summary.sort(key=lambda x: (-x["mean_preference"], -x["mean_male"]))
    report = {
        "kind": "posthoc_exploratory_parameter_search",
        "not_pre_registered": True,
        "baseline_unchanged": True,
        "reversal_mV": REVERSAL,
        "seeds": list(SEEDS),
        "seconds": time.perf_counter() - started,
        "results": rows,
        "summary": summary,
        "flipped_conditions": [s for s in summary if s["all_seeds_male_gt_female"]],
    }
    out = Path("exploratory-male-bias.json")
    out.write_text(json.dumps(report, indent=2))
    print("=== SUMMARY (best preference first) ===", flush=True)
    for s in summary:
        print(
            json.dumps(
                {
                    "condition": s["condition"],
                    "mean_M": s["mean_male"],
                    "mean_F": s["mean_female"],
                    "pref": s["mean_preference"],
                    "wins": s["seeds_male_gt_female"],
                }
            ),
            flush=True,
        )
    print("COMPLETE", report["seconds"], flush=True)


if __name__ == "__main__":
    main()
