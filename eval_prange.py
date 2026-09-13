"""Evaluate whether membrane integration in `_simulate` can use prange.

Does NOT modify the production kernel. Implements an experimental variant that
parallelizes ONLY the per-neuron subthreshold update (race-free), and keeps the
synaptic scatter serial (shared-write race on pending_e/pending_i).

Reports bitwise equality vs production `_simulate` and wall times.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from numba import njit, prange, set_num_threads

from orientation.intervention import none, output_silence
from orientation.simulate import (
    BIN_MS,
    DELAY_MS,
    DT_MS,
    EXTERNAL_MV,
    MEMBRANE_TAU_MS,
    REFRACTORY_MS,
    REST_MV,
    SYN_TAU_MS,
    THRESH_MV,
    WEIGHT_SCALE,
    _simulate,
    apply_sign_overrides,
    load_prepared,
    sensory_event_schedule,
)


@njit(cache=True, parallel=True)
def _simulate_prange_membrane(
    indptr,
    dest,
    weights,
    signs,
    sensory_inputs,
    sensory_events,
    output_mask,
    output_gain,
    tonic,
    activation_hz,
    activation_target,
    ei,
    dt,
):
    """Same equations as `_simulate`; membrane loop uses prange; scatter serial."""
    n = signs.shape[0]
    steps = sensory_events.shape[0]
    delay = int(round(DELAY_MS / dt))
    ref = int(round(REFRACTORY_MS / dt))
    v = np.full(n, REST_MV)
    ge = np.zeros(n)
    hi = np.zeros(n)
    release = np.zeros(n, np.int32)
    pending_e = np.zeros((delay + 1, n))
    pending_i = np.zeros((delay + 1, n))
    stim_mask = np.zeros(n, np.bool_)
    act_mask = np.zeros(n, np.bool_)
    for j in range(sensory_inputs.shape[0]):
        stim_mask[sensory_inputs[j]] = True
    for j in range(activation_target.shape[0]):
        act_mask[activation_target[j]] = True
    n_bins = int(np.ceil(steps * dt / BIN_MS))
    counts = np.zeros((n_bins, n), np.int32)
    eg = np.exp(-dt / SYN_TAU_MS)
    vmin = REST_MV
    vmax = REST_MV
    act_p = activation_hz * dt / 1000.0
    local_vmin = np.full(1, REST_MV)
    local_vmax = np.full(1, REST_MV)
    for step in range(steps):
        t = step * dt
        slot = step % (delay + 1)
        # Race-free: each i reads/writes only its own columns of pending/v/ge/hi.
        for i in prange(n):
            inc_e = pending_e[slot, i]
            inc_i = pending_i[slot, i]
            pending_e[slot, i] = 0.0
            pending_i[slot, i] = 0.0
            if step >= release[i]:
                ge[i] += inc_e
                hi[i] += inc_i
                rate = (1.0 + hi[i]) / MEMBRANE_TAU_MS
                decay = np.exp(-rate * dt)
                drive = tonic[i] if act_mask[i] else 0.0
                equilibrium = (REST_MV + hi[i] * ei + drive) / (1.0 + hi[i])
                denom = rate - 1.0 / SYN_TAU_MS
                if abs(denom) > 1e-10:
                    coupling = (eg - decay) / (MEMBRANE_TAU_MS * denom)
                else:
                    coupling = dt * decay / MEMBRANE_TAU_MS
                v[i] = equilibrium + (v[i] - equilibrium) * decay + ge[i] * coupling
                ge[i] *= eg
                hi[i] *= eg
            # reduction over v for extrema would race if written to shared scalars;
            # recompute serially below from v (cheap relative to scatter).
        for j in range(sensory_inputs.shape[0]):
            if sensory_events[step, j]:
                v[sensory_inputs[j]] += EXTERNAL_MV
        if act_p > 0.0:
            for j in range(activation_target.shape[0]):
                if np.random.random() < act_p:
                    v[activation_target[j]] += EXTERNAL_MV
        # Serial spike detection + scatter: scatter has shared-write races.
        for i in range(n):
            if v[i] < vmin:
                vmin = v[i]
            if v[i] > vmax:
                vmax = v[i]
            if step >= release[i] and v[i] > THRESH_MV:
                bin_i = int(t / BIN_MS)
                if bin_i >= n_bins:
                    bin_i = n_bins - 1
                counts[bin_i, i] += 1
                v[i] = REST_MV
                ge[i] = 0.0
                hi[i] = 0.0
                release[i] = step + (0 if stim_mask[i] else ref)
                if signs[i] == 0:
                    continue
                gain = output_gain if output_mask[i] else 1.0
                if gain == 0.0:
                    continue
                slot2 = (step + delay) % (delay + 1)
                for e in range(indptr[i], indptr[i + 1]):
                    target = dest[e]
                    w = weights[e] * WEIGHT_SCALE * gain
                    if signs[i] > 0:
                        pending_e[slot2, target] += w
                    else:
                        pending_i[slot2, target] += w / (REST_MV - ei)
    return counts, vmin, vmax


def pack(net, interv, input_name, seed):
    signs = apply_sign_overrides(net)
    inputs, events = sensory_event_schedule(net, input_name, seed)
    packed = interv.packed(net.n)
    act_target = (
        interv.target if interv.kind == "activation" else np.zeros(0, np.int32)
    )
    return (
        net.graph.indptr,
        net.graph.indices,
        net.graph.data,
        signs,
        inputs if input_name != "no_input" else np.zeros(0, np.int32),
        events,
        packed["output_mask"],
        float(packed["output_gain"]),
        packed["tonic"],
        float(packed["activation_hz"]),
        act_target,
        -70.0,
        DT_MS,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "prange-eval.json",
    )
    args = parser.parse_args()

    net = load_prepared(args.data)
    interv = output_silence(net.groups["mAL"], name="mAL_output_silence")
    args_pack = pack(net, interv, "candidate_male", 11)

    # Warmup JIT
    _simulate(*args_pack)
    _simulate_prange_membrane(*args_pack)

    results = []
    for threads in (1, 2, 4, 8):
        set_num_threads(threads)
        # equality
        c0, v0, w0 = _simulate(*args_pack)
        c1, v1, w1 = _simulate_prange_membrane(*args_pack)
        counts_equal = bool(np.array_equal(c0, c1))
        extrema_equal = bool(v0 == v1 and w0 == w1)
        equal = counts_equal and extrema_equal
        # timing
        t_prod = []
        t_par = []
        for _ in range(args.repeats):
            t = time.perf_counter()
            _simulate(*args_pack)
            t_prod.append(time.perf_counter() - t)
            t = time.perf_counter()
            _simulate_prange_membrane(*args_pack)
            t_par.append(time.perf_counter() - t)
        row = {
            "threads": threads,
            "counts_equal": counts_equal,
            "extrema_equal": extrema_equal,
            "bitwise_equal": equal,
            "vmin_prod": v0,
            "vmin_prange": v1,
            "vmax_prod": w0,
            "vmax_prange": w1,
            "P1_spikes_prod": int(c0[5:, net.groups["P1_readout"]].sum()),
            "P1_spikes_prange": int(c1[5:, net.groups["P1_readout"]].sum()),
            "network_spikes_prod": int(c0.sum()),
            "network_spikes_prange": int(c1.sum()),
            "prod_median_s": float(np.median(t_prod)),
            "prange_median_s": float(np.median(t_par)),
            "speedup": float(np.median(t_prod) / np.median(t_par)),
        }
        results.append(row)
        print(json.dumps(row), flush=True)

    report = {
        "note": (
            "Experimental probe only. Production `_simulate` unchanged. "
            "prange applies only to subthreshold membrane update; synaptic "
            "scatter remains serial (shared-write race)."
        ),
        "results": results,
        "recommendation": None,
    }
    if all(r["counts_equal"] for r in results):
        best = max(results, key=lambda r: r["speedup"])
        if best["speedup"] > 1.05:
            extrema_ok = all(r["extrema_equal"] for r in results)
            report["recommendation"] = (
                "Spike counts match production exactly. "
                + (
                    "Extrema also match."
                    if extrema_ok
                    else "Voltage extrema differ (scan point); diagnostics only."
                )
                + f" Measured up to {best['speedup']:.2f}x at threads={best['threads']} "
                "if membrane-only prange is adopted; keep synaptic scatter serial."
            )
        else:
            report["recommendation"] = (
                "prange membrane update is count-identical but no reliable speedup; "
                "keep serial kernel."
            )
    else:
        report["recommendation"] = (
            "Do NOT adopt prange variant: spike counts differ from production kernel."
        )
    args.out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
