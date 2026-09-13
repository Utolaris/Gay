"""Experimental kernels for E3–E8. Production simulate.py is unchanged."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pyarrow.feather as feather
from numba import njit
from scipy.sparse import csr_matrix, load_npz

from .simulate import (
    BIN_MS,
    DELAY_MS,
    DT_MS,
    EXTERNAL_MV,
    MEMBRANE_TAU_MS,
    P1_BODY_IDS,
    REFRACTORY_MS,
    REST_MV,
    SENSORY_TOTAL_HZ,
    STIM_END_MS,
    STIM_START_MS,
    SYN_TAU_MS,
    THRESH_MV,
    VAB3_BODY_IDS,
    WEIGHT_SCALE,
    DURATION_MS,
    Network,
)


@njit(cache=True)
def simulate_ex(
    indptr,
    dest,
    weights,
    signs,
    sensory_inputs,
    sensory_events,
    neuron_gain,
    tonic,
    ei,
    use_arousal,
    arousal_global,
    arousal_s_max,
    arousal_alpha,
    arousal_tau_ms,
    male_relay_mask,
):
    """Bounded LIF + per-neuron gain/tonic + optional arousal state.

    Arousal S multiplies outgoing weights: globally, or only when presynaptic
    male_relay_mask[i] is true. Spike scatter remains serial (shared writes).
    """
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
    vmin = REST_MV
    vmax = REST_MV
    S = 1.0
    a_decay = np.exp(-DT_MS / arousal_tau_ms)
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
            if v[i] < vmin:
                vmin = v[i]
            if v[i] > vmax:
                vmax = v[i]
        for j in range(sensory_inputs.shape[0]):
            if sensory_events[step, j]:
                v[sensory_inputs[j]] += EXTERNAL_MV
        n_spikes_step = 0
        for i in range(n):
            if step >= release[i] and v[i] > THRESH_MV:
                n_spikes_step += 1
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
                g = neuron_gain[i]
                if use_arousal:
                    if arousal_global or male_relay_mask[i]:
                        g = g * S
                slot2 = (step + delay) % (delay + 1)
                for e in range(indptr[i], indptr[i + 1]):
                    target = dest[e]
                    w = weights[e] * WEIGHT_SCALE * g
                    if signs[i] > 0:
                        pending_e[slot2, target] += w
                    else:
                        pending_i[slot2, target] += w / (REST_MV - ei)
        if use_arousal and n_spikes_step > 0:
            S = S + arousal_alpha * n_spikes_step
            if S > arousal_s_max:
                S = arousal_s_max
        if use_arousal:
            S = 1.0 + (S - 1.0) * a_decay
    return counts, vmin, vmax


def load_ex_network(data_dir: Path, e8_mult: float = 1.0) -> Network:
    root = Path(data_dir)
    meta = json.loads((root / "full-meta.json").read_text())
    nodes = np.load(root / "full-nodes.npz")
    ids = np.asarray(nodes["ids"])
    signs = np.asarray(nodes["signs"], dtype=np.int8).copy()
    graph = load_npz(root / "full-graph.npz").tocsr()
    p1 = np.searchsorted(ids, P1_BODY_IDS)
    ef = np.searchsorted(ids, VAB3_BODY_IDS)
    if e8_mult != 1.0:
        graph = graph.copy()
        p1set = set(int(x) for x in p1)
        indptr, indices, data = graph.indptr, graph.indices, graph.data.copy()
        for i in ef:
            for e in range(indptr[i], indptr[i + 1]):
                if int(indices[e]) in p1set:
                    data[e] = data[e] * e8_mult
        graph = csr_matrix((data, indices, indptr), shape=graph.shape)
    groups = {k: np.asarray(v, dtype=np.int32) for k, v in meta["groups"].items()}
    groups["vAB3"] = ef.astype(np.int32)
    groups["P1_readout"] = p1.astype(np.int32)
    signs[groups["vAB3"]] = 1
    signs[groups["candidate_female"]] = 1
    return Network(graph=graph, signs=signs, ids=ids, groups=groups)


def type_indices(net: Network, data_dir: Path, types: list[str]) -> np.ndarray:
    ann = feather.read_table(Path(data_dir) / "annotations.feather").to_pylist()
    by_id = {r["bodyId"]: r for r in ann}
    want = set(types)
    return np.asarray(
        [
            i
            for i, bid in enumerate(net.ids)
            if (by_id.get(int(bid)) or {}).get("type") in want
        ],
        dtype=np.int32,
    )


def make_events_with_maps(
    input_name: str,
    seed: int,
    encoding: str,
    full_idx: np.ndarray,
    type_maps: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    steps = int(round(DURATION_MS / DT_MS))
    s0 = int(round(STIM_START_MS / DT_MS))
    s1 = int(round(STIM_END_MS / DT_MS))
    if encoding == "full":
        inputs = full_idx
    elif encoding == "half_random":
        inputs = np.sort(full_idx)[: max(1, len(full_idx) // 2)]
    elif encoding == "type_only_a":
        key = "male_a" if input_name == "candidate_male" else "female_a"
        inputs = type_maps[key]
    elif encoding == "type_only_b":
        key = "male_b" if input_name == "candidate_male" else "female_b"
        inputs = type_maps[key]
    else:
        raise ValueError(encoding)
    inputs = np.asarray(inputs, dtype=np.int32)
    rate = SENSORY_TOTAL_HZ / max(len(inputs), 1)
    ev = np.zeros((steps, len(inputs)), dtype=np.bool_)
    ev[s0:s1] = np.random.RandomState(seed).random_sample((s1 - s0, len(inputs))) < (
        rate * DT_MS / 1000.0
    )
    return inputs, ev


def run_ex(
    net: Network,
    input_name: str,
    seed: int,
    *,
    gain: np.ndarray | None = None,
    tonic: np.ndarray | None = None,
    ei: float = -70.0,
    use_arousal: bool = False,
    arousal_global: bool = True,
    arousal_s_max: float = 2.0,
    arousal_alpha: float = 0.02,
    arousal_tau_ms: float = 50.0,
    male_relay_mask: np.ndarray | None = None,
    inputs: np.ndarray | None = None,
    events: np.ndarray | None = None,
    encoding: str = "full",
    type_maps: dict[str, np.ndarray] | None = None,
) -> dict:
    n = net.n
    if inputs is None or events is None:
        if encoding == "full" and input_name != "no_input":
            from .simulate import sensory_event_schedule

            inputs, events = sensory_event_schedule(net, input_name, seed)
        elif input_name == "no_input":
            steps = int(round(DURATION_MS / DT_MS))
            inputs = np.zeros(0, np.int32)
            events = np.zeros((steps, 0), dtype=np.bool_)
        else:
            if type_maps is None:
                raise ValueError("type_maps required for non-full encoding")
            inputs, events = make_events_with_maps(
                input_name, seed, encoding, net.groups[input_name], type_maps
            )
    g = np.ones(n, dtype=np.float64) if gain is None else gain
    ton = np.zeros(n, dtype=np.float64) if tonic is None else tonic
    relay = np.zeros(n, dtype=np.bool_) if male_relay_mask is None else male_relay_mask
    counts, vmin, vmax = simulate_ex(
        net.graph.indptr,
        net.graph.indices,
        net.graph.data,
        net.signs,
        inputs,
        events,
        g,
        ton,
        float(ei),
        use_arousal,
        arousal_global,
        float(arousal_s_max),
        float(arousal_alpha),
        float(arousal_tau_ms),
        relay,
    )
    if vmin < ei - 1e-8:
        raise AssertionError(f"voltage {vmin} below reversal {ei}")
    p1 = net.groups["P1_readout"]
    broad = net.groups["P1_related"]
    window = counts[5:]
    return {
        "input": input_name,
        "seed": int(seed),
        "P1_spikes": int(window[:, p1].sum()),
        "broad_P1_spikes": int(window[:, broad].sum()),
        "mAL_spikes": int(window[:, net.groups["mAL"]].sum()),
        "network_spikes": int(counts.sum()),
        "active_neurons": int(np.any(counts, axis=0).sum()),
        "events": int(events.sum()),
        "voltage_min_mV": float(vmin),
        "voltage_max_mV": float(vmax),
        "n_input_cells": int(len(inputs)),
    }
