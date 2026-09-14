"""P1-targeted synaptic current tracking (delivered weight flux).

Does not modify production simulate.py. Same frozen LIF constants.
Attribution uses actual arrivals onto P1, not source spike counts.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyarrow.feather as feather
from numba import njit

from .simulate import (
    BIN_MS,
    DELAY_MS,
    DT_MS,
    EXTERNAL_MV,
    MEMBRANE_TAU_MS,
    REFRACTORY_MS,
    REST_MV,
    STIM_END_MS,
    STIM_START_MS,
    SYN_TAU_MS,
    THRESH_MV,
    DURATION_MS,
    WEIGHT_SCALE,
    Network,
    sensory_event_schedule,
)


@dataclass
class P1Map:
    p1: np.ndarray  # graph indices of 8 P1 cells
    p1_local: np.ndarray  # graph index -> 0..7 or -1
    type_ids: np.ndarray  # per neuron, type id or -1
    type_names: list[str]  # id -> type name
    is_inh_type: np.ndarray  # per type id, whether sources are inhibitory (sign -1)
    # any neuron with sign==0 never delivers; types still named


def build_p1_map(net: Network, data_dir: Path) -> P1Map:
    ann = feather.read_table(Path(data_dir) / "annotations.feather").to_pylist()
    by_id = {r["bodyId"]: r for r in ann}
    p1 = np.asarray(net.groups["P1_readout"], dtype=np.int32)
    p1_local = np.full(net.n, -1, dtype=np.int32)
    p1_local[p1] = np.arange(p1.size, dtype=np.int32)
    p1set = set(int(x) for x in p1)

    # types of all presynaptic partners of P1 + the P1 cells themselves
    indptr, indices = net.graph.indptr, net.graph.indices
    # walk incoming via CSR transpose = iterate all edges to p1
    partner_types: dict[str, None] = {}
    # Use CSC-equivalent: scan all rows for edges into p1 (expensive) —
    # better: build from CSR by checking each source's outs. Only need types
    # that actually deliver; we assign type_id to every neuron anyway.
    type_names: list[str] = ["__none__"]
    type_index: dict[str, int] = {"__none__": 0}
    type_ids = np.zeros(net.n, dtype=np.int32)
    for i, bid in enumerate(net.ids):
        t = (by_id.get(int(bid)) or {}).get("type") or f"body_{int(bid)}"
        if t not in type_index:
            type_index[t] = len(type_names)
            type_names.append(t)
        type_ids[i] = type_index[t]
    is_inh = np.zeros(len(type_names), dtype=np.bool_)
    # type-level sign from first member
    seen = np.zeros(len(type_names), dtype=np.bool_)
    for i in range(net.n):
        tid = type_ids[i]
        if not seen[tid]:
            seen[tid] = True
            is_inh[tid] = net.signs[i] < 0
    return P1Map(
        p1=p1,
        p1_local=p1_local,
        type_ids=type_ids,
        type_names=type_names,
        is_inh_type=is_inh,
    )


@njit(cache=True)
def _simulate_p1_currents(
    indptr,
    dest,
    weights,
    signs,
    sensory_inputs,
    sensory_events,
    output_gain,
    tonic,
    ei,
    p1_cells,
    p1_local,
    type_ids,
    n_types,
    n_bins,
):
    """Return P1 traces and per-type delivered flux.

    Arrays:
      vm[bin, 8], ge[bin, 8], hi[bin, 8]
      e_flux[bin, 8], i_flux[bin, 8]
      e_by_type[bin, n_types], i_by_type[bin, n_types]  # summed over P1 cells
      first_spike_bin[8]  # -1 if none
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
    eg = np.exp(-DT_MS / SYN_TAU_MS)
    act_p = 0.0

    np1 = p1_cells.shape[0]
    vm = np.zeros((n_bins, np1))
    ge_tr = np.zeros((n_bins, np1))
    hi_tr = np.zeros((n_bins, np1))
    e_flux = np.zeros((n_bins, np1))
    i_flux = np.zeros((n_bins, np1))
    e_by_type = np.zeros((n_bins, n_types))
    i_by_type = np.zeros((n_bins, n_types))
    first_spike_bin = np.full(np1, -1, np.int32)
    p1_counts = np.zeros(np1, np.int32)

    # per-step flux accumulators then fold into bin
    step_e = np.zeros(np1)
    step_i = np.zeros(np1)
    # type flux per step is expensive; accumulate directly into current bin

    for step in range(steps):
        t = step * DT_MS
        bin_i = int(t / BIN_MS)
        if bin_i >= n_bins:
            bin_i = n_bins - 1
        slot = step % (delay + 1)
        # arrivals: add to ge/hi and count delivered flux for P1
        for i in range(n):
            inc_e = pending_e[slot, i]
            inc_i = pending_i[slot, i]
            pending_e[slot, i] = 0.0
            pending_i[slot, i] = 0.0
            if inc_e != 0.0 or inc_i != 0.0:
                loc = p1_local[i]
                if loc >= 0:
                    if inc_e != 0.0:
                        e_flux[bin_i, loc] += inc_e
                    if inc_i != 0.0:
                        i_flux[bin_i, loc] += inc_i
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
        # record Vm for P1 at end of voltage update (before possible spike reset)
        for li in range(np1):
            gi = p1_cells[li]
            vm[bin_i, li] = v[gi]
            ge_tr[bin_i, li] = ge[gi]
            hi_tr[bin_i, li] = hi[gi]
        # spikes + scatter
        for i in range(n):
            if step >= release[i] and v[i] > THRESH_MV:
                loc = p1_local[i]
                if loc >= 0:
                    p1_counts[loc] += 1
                    if first_spike_bin[loc] < 0:
                        first_spike_bin[loc] = bin_i
                v[i] = REST_MV
                ge[i] = 0.0
                hi[i] = 0.0
                release[i] = step + (0 if stim_mask[i] else ref)
                if signs[i] == 0:
                    continue
                gain = output_gain[i]
                if gain == 0.0:
                    continue
                sg = signs[i]
                tid = type_ids[i]
                slot2 = (step + delay) % (delay + 1)
                # bin for delivery is (step+delay)
                tdel = t + DELAY_MS
                bdel = int(tdel / BIN_MS)
                if bdel >= n_bins:
                    bdel = n_bins - 1
                for e in range(indptr[i], indptr[i + 1]):
                    target = dest[e]
                    w = weights[e] * WEIGHT_SCALE * gain
                    loct = p1_local[target]
                    if sg > 0:
                        pending_e[slot2, target] += w
                        if loct >= 0:
                            e_by_type[bdel, tid] += w
                    else:
                        pending_i[slot2, target] += w / (REST_MV - ei)
                        if loct >= 0:
                            # store unscaled anatomical*gain weight for ranking
                            # (same units as e_by_type); actual conductance uses
                            # the / (REST-ei) factor in the kernel.
                            i_by_type[bdel, tid] += w
    return vm, ge_tr, hi_tr, e_flux, i_flux, e_by_type, i_by_type, first_spike_bin, p1_counts


def run_p1_current_trial(
    net: Network,
    pmap: P1Map,
    signs: np.ndarray,
    input_name: str,
    seed: int,
    *,
    gain: np.ndarray | None = None,
    tonic: np.ndarray | None = None,
) -> dict:
    n = net.n
    steps = int(round(DURATION_MS / DT_MS))
    n_bins = int(np.ceil(steps * DT_MS / BIN_MS))
    if input_name == "no_input":
        inputs = np.zeros(0, np.int32)
        events = np.zeros((steps, 0), np.bool_)
    else:
        inputs, events = sensory_event_schedule(net, input_name, seed)
    g = np.ones(n, dtype=np.float64) if gain is None else gain
    ton = np.zeros(n, dtype=np.float64) if tonic is None else tonic
    vm, ge, hi, e_flux, i_flux, e_by_type, i_by_type, first_spike, p1_counts = (
        _simulate_p1_currents(
            net.graph.indptr,
            net.graph.indices,
            net.graph.data,
            signs,
            inputs,
            events,
            g,
            ton,
            -70.0,
            pmap.p1,
            pmap.p1_local,
            pmap.type_ids,
            len(pmap.type_names),
            n_bins,
        )
    )
    return {
        "input": input_name,
        "seed": int(seed),
        "bin_ms": BIN_MS,
        "n_bins": n_bins,
        "vm": vm.tolist(),
        "ge": ge.tolist(),
        "hi": hi.tolist(),
        "e_flux": e_flux.tolist(),
        "i_flux": i_flux.tolist(),
        "e_by_type": e_by_type.tolist(),
        "i_by_type": i_by_type.tolist(),
        "first_spike_bin": first_spike.tolist(),
        "p1_cell_spikes": p1_counts.tolist(),
        "p1_spikes": int(p1_counts.sum()),
        "n_p1_active": int((p1_counts > 0).sum()),
    }
