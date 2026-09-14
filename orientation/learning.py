"""DA-gated KC→MBON plasticity + three-phase reward learning.

Does not modify production simulate.py. Frozen LIF constants imported.
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

# Pre-registered (REWARD_PROTOCOL.md)
ETA = 0.05
TAU_ELIG_MS = 50.0
W_FLOOR = 0.25
PAM_HZ = 20.0
KC_ACCESS_HZ = 20.0
KC_ACCESS_N = 200
N_EPOCHS = 20

KC_TYPES = (
    "KCg-m",
    "KCab-s",
    "KCab-m",
    "KCab-c",
    "KCa'b'-ap2",
    "KCg-d",
    "KCa'b'-m",
    "KCa'b'-ap1",
    "KCab-p",
)

MBON_TYPES = tuple([f"MBON{i:02d}" for i in range(1, 36)]) + (
    "MBON15-like",
    "MBON25-like",
    "MBON17-like",
)


@dataclass
class PlasticNet:
    net: Network
    kc: np.ndarray
    mbon: np.ndarray
    pam: np.ndarray
    plast_indptr: np.ndarray
    plast_dest: np.ndarray
    plast_w0: np.ndarray
    plast_w: np.ndarray
    kc_pos: np.ndarray


def type_indices(net: Network, data_dir: Path, types) -> np.ndarray:
    ann = feather.read_table(Path(data_dir) / "annotations.feather").to_pylist()
    by_id = {r["bodyId"]: r for r in ann}
    want = set(types)
    return np.asarray(
        [i for i, bid in enumerate(net.ids) if (by_id.get(int(bid)) or {}).get("type") in want],
        dtype=np.int32,
    )


def build_plastic_net(net: Network, data_dir: Path) -> PlasticNet:
    kc = type_indices(net, data_dir, KC_TYPES)
    mbon = type_indices(net, data_dir, MBON_TYPES)
    pam = type_indices(net, data_dir, [f"PAM{i:02d}" for i in range(1, 16)])
    mbon_set = set(int(x) for x in mbon)
    kc_pos = np.full(net.n, -1, dtype=np.int32)
    kc_pos[kc] = np.arange(kc.size, dtype=np.int32)

    indptr, indices, data = net.graph.indptr, net.graph.indices, net.graph.data
    dest_list: list[int] = []
    w0_list: list[float] = []
    indptr_l = [0]
    for gi in kc:
        s, e = indptr[gi], indptr[gi + 1]
        for t, w in zip(indices[s:e], data[s:e]):
            t = int(t)
            if t in mbon_set:
                dest_list.append(t)
                w0_list.append(float(w))
        indptr_l.append(len(dest_list))
    return PlasticNet(
        net=net,
        kc=kc,
        mbon=mbon,
        pam=pam,
        plast_indptr=np.asarray(indptr_l, dtype=np.int32),
        plast_dest=np.asarray(dest_list, dtype=np.int32),
        plast_w0=np.asarray(w0_list, dtype=np.float64),
        plast_w=np.asarray(w0_list, dtype=np.float64).copy(),
        kc_pos=kc_pos,
    )


@njit(cache=True)
def _simulate_learn(
    indptr,
    dest,
    weights,
    signs,
    sensory_inputs,
    sensory_events,
    output_gain,
    tonic,
    ei,
    kc_pos,
    plast_indptr,
    plast_dest,
    plast_w,
    plast_w0,
    elig,
    da_flag,
    learn_on,
    act_target,
    act_hz,
    eta,
    tau_elig,
    w_floor,
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
    act_mask = np.zeros(n, np.bool_)
    for j in range(sensory_inputs.shape[0]):
        stim_mask[sensory_inputs[j]] = True
    for j in range(act_target.shape[0]):
        act_mask[act_target[j]] = True
    n_bins = int(np.ceil(steps * DT_MS / BIN_MS))
    counts = np.zeros((n_bins, n), np.int32)
    eg = np.exp(-DT_MS / SYN_TAU_MS)
    e_decay = np.exp(-DT_MS / tau_elig)
    vmin = REST_MV
    vmax = REST_MV
    act_p = act_hz * DT_MS / 1000.0

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
                drive = tonic[i] if act_mask[i] else 0.0
                equilibrium = (REST_MV + hi[i] * ei + drive) / (1.0 + hi[i])
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
        if act_p > 0.0:
            for j in range(act_target.shape[0]):
                if np.random.random() < act_p:
                    v[act_target[j]] += EXTERNAL_MV
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
                k = kc_pos[i]
                if k >= 0:
                    elig[k] = 1.0
                if signs[i] == 0:
                    continue
                gain = output_gain[i]
                if gain == 0.0:
                    continue
                slot2 = (step + delay) % (delay + 1)
                s0 = indptr[i]
                s1 = indptr[i + 1]
                if k >= 0:
                    p0 = plast_indptr[k]
                    p1 = plast_indptr[k + 1]
                    for e in range(s0, s1):
                        target = dest[e]
                        w = weights[e] * WEIGHT_SCALE * gain
                        for pe in range(p0, p1):
                            if plast_dest[pe] == target:
                                w = plast_w[pe] * WEIGHT_SCALE * gain
                                break
                        if signs[i] > 0:
                            pending_e[slot2, target] += w
                        else:
                            pending_i[slot2, target] += w / (REST_MV - ei)
                else:
                    for e in range(s0, s1):
                        target = dest[e]
                        w = weights[e] * WEIGHT_SCALE * gain
                        if signs[i] > 0:
                            pending_e[slot2, target] += w
                        else:
                            pending_i[slot2, target] += w / (REST_MV - ei)
        if learn_on and da_flag[step]:
            for k in range(elig.shape[0]):
                e0 = elig[k]
                if e0 <= 0.0:
                    continue
                for pe in range(plast_indptr[k], plast_indptr[k + 1]):
                    floor = w_floor * plast_w0[pe]
                    nw = plast_w[pe] - eta * e0
                    if nw < floor:
                        nw = floor
                    if nw > plast_w0[pe]:
                        nw = plast_w0[pe]
                    plast_w[pe] = nw
        for k in range(elig.shape[0]):
            elig[k] *= e_decay
    return counts, vmin, vmax


def make_da_flag(n_steps: int, on: bool) -> np.ndarray:
    da = np.zeros(n_steps, np.bool_)
    if on:
        da[int(STIM_START_MS / DT_MS) : int(STIM_END_MS / DT_MS)] = True
    return da


def run_learn_trial(
    pnet: PlasticNet,
    signs: np.ndarray,
    input_name: str,
    seed: int,
    *,
    da_on: bool = False,
    learn_on: bool = True,
    gain: np.ndarray | None = None,
    tonic: np.ndarray | None = None,
    kc_access: bool = False,
    reward_hz: float = 0.0,
    events: np.ndarray | None = None,
    inputs: np.ndarray | None = None,
) -> tuple[np.ndarray, dict]:
    net = pnet.net
    n = net.n
    steps = int(round(DURATION_MS / DT_MS))
    if inputs is None or events is None:
        if input_name == "no_input":
            inputs = np.zeros(0, np.int32)
            events = np.zeros((steps, 0), np.bool_)
        else:
            inputs, events = sensory_event_schedule(net, input_name, seed)
    g = np.ones(n, dtype=np.float64) if gain is None else gain
    ton = np.zeros(n, dtype=np.float64) if tonic is None else tonic

    parts: list[np.ndarray] = []
    act_hz = 0.0
    if reward_hz > 0.0:
        parts.append(pnet.pam)
        act_hz = float(reward_hz)
    if kc_access:
        rng = np.random.RandomState(seed + 77)
        sub = np.sort(
            rng.choice(pnet.kc, size=min(KC_ACCESS_N, pnet.kc.size), replace=False)
        ).astype(np.int32)
        parts.append(sub)
        if act_hz == 0.0:
            act_hz = KC_ACCESS_HZ
    if parts:
        act_target = np.unique(np.concatenate(parts)).astype(np.int32)
    else:
        act_target = np.zeros(0, dtype=np.int32)

    da_flag = make_da_flag(steps, da_on)
    elig = np.zeros(pnet.kc.size, dtype=np.float64)
    counts, vmin, vmax = _simulate_learn(
        net.graph.indptr,
        net.graph.indices,
        net.graph.data,
        signs,
        inputs,
        events,
        g,
        ton,
        -70.0,
        pnet.kc_pos,
        pnet.plast_indptr,
        pnet.plast_dest,
        pnet.plast_w,
        pnet.plast_w0,
        elig,
        da_flag,
        learn_on,
        act_target,
        act_hz,
        ETA,
        TAU_ELIG_MS,
        W_FLOOR,
    )
    if vmin < -70.0 - 1e-8:
        raise AssertionError(f"voltage {vmin} below reversal")
    return counts, {
        "events": int(events.sum()) if events.size else 0,
        "event_seed": int(seed),
        "vmin": float(vmin),
        "vmax": float(vmax),
        "act_n": int(act_target.size),
    }


def reset_weights(pnet: PlasticNet) -> None:
    pnet.plast_w[:] = pnet.plast_w0


def weight_stats(pnet: PlasticNet) -> dict:
    w, w0 = pnet.plast_w, pnet.plast_w0
    if w.size == 0:
        return {"n_edges": 0, "mean_ratio": 0.0, "frac_depressed": 0.0, "frac_at_floor": 0.0, "total_w": 0.0, "total_w0": 0.0}
    ratio = w / np.maximum(w0, 1e-12)
    return {
        "n_edges": int(w.size),
        "mean_ratio": float(ratio.mean()),
        "frac_depressed": float((ratio < 0.999).mean()),
        "frac_at_floor": float((ratio <= W_FLOOR + 1e-6).mean()),
        "total_w": float(w.sum()),
        "total_w0": float(w0.sum()),
    }


def class_spikes(counts: np.ndarray, idx: np.ndarray) -> int:
    return int(counts[5:, idx].sum())
