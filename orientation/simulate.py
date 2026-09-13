"""Bounded-conductance LIF kernel with class-level intervention hooks.

Equations follow experiment/followup/bounded_routes.py so baseline WT/mAL-silence
results are comparable to the upstream follow-up without importing or editing it.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass

import numpy as np
from numba import njit
from scipy.sparse import csr_matrix

from .intervention import Intervention, none

# Frozen LIF / protocol constants (do not retune).
DT_MS = 0.1
DELAY_MS = 1.8
REFRACTORY_MS = 2.2
REST_MV = -52.0
THRESH_MV = -45.0
MEMBRANE_TAU_MS = 20.0
SYN_TAU_MS = 5.0
WEIGHT_SCALE = 0.275
EXTERNAL_MV = 68.75
DURATION_MS = 300.0
STIM_START_MS = 50.0
STIM_END_MS = 250.0
BIN_MS = 10.0
DEFAULT_REVERSAL_MV = -70.0
SENSORY_TOTAL_HZ = 5550.0
P1_BODY_IDS = np.array(
    [12442, 16719, 17867, 20117, 20803, 23968, 519518, 522419], dtype=np.int64
)
VAB3_BODY_IDS = np.array([11998, 13341, 13693, 512498], dtype=np.int64)


@njit(cache=True)
def _simulate(
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
    for step in range(steps):
        t = step * dt
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
            if v[i] < vmin:
                vmin = v[i]
            if v[i] > vmax:
                vmax = v[i]
        for j in range(sensory_inputs.shape[0]):
            if sensory_events[step, j]:
                v[sensory_inputs[j]] += EXTERNAL_MV
        if act_p > 0.0:
            for j in range(activation_target.shape[0]):
                if np.random.random() < act_p:
                    v[activation_target[j]] += EXTERNAL_MV
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


@dataclass(frozen=True)
class Network:
    graph: csr_matrix
    signs: np.ndarray
    ids: np.ndarray
    groups: dict[str, np.ndarray]

    @property
    def n(self) -> int:
        return self.graph.shape[0]


@dataclass(frozen=True)
class TrialResult:
    counts: np.ndarray
    voltage_min_mV: float
    voltage_max_mV: float
    events: int
    event_sha256: str
    input_name: str
    seed: int
    intervention: str
    reversal_mV: float


def load_prepared(data_dir) -> Network:
    from pathlib import Path

    import json

    from scipy.sparse import load_npz

    root = Path(data_dir)
    meta = json.loads((root / "full-meta.json").read_text())
    nodes = np.load(root / "full-nodes.npz")
    graph = load_npz(root / "full-graph.npz").tocsr()
    groups = {k: np.asarray(v, dtype=np.int32) for k, v in meta["groups"].items()}
    groups["vAB3"] = np.searchsorted(nodes["ids"], VAB3_BODY_IDS).astype(np.int32)
    groups["P1_readout"] = np.searchsorted(nodes["ids"], P1_BODY_IDS).astype(np.int32)
    if not np.array_equal(nodes["ids"][groups["P1_readout"]], P1_BODY_IDS):
        raise ValueError("P1 body IDs not found in prepared graph")
    return Network(
        graph=graph,
        signs=np.asarray(nodes["signs"], dtype=np.int8),
        ids=np.asarray(nodes["ids"]),
        groups=groups,
    )


def apply_sign_overrides(network: Network) -> np.ndarray:
    """Upstream bounded sensory sign map, identical across inputs."""
    signs = network.signs.copy()
    signs[network.groups["vAB3"]] = 1
    signs[network.groups["candidate_female"]] = 1
    return signs


def sensory_event_schedule(
    network: Network, input_name: str, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    """Generate one matched Bernoulli schedule for male/female/no-input."""
    steps = int(round(DURATION_MS / DT_MS))
    stim_start = int(round(STIM_START_MS / DT_MS))
    stim_end = int(round(STIM_END_MS / DT_MS))
    if input_name == "no_input":
        return (
            np.zeros(0, dtype=np.int32),
            np.zeros((steps, 0), dtype=np.bool_),
        )
    if input_name not in ("candidate_male", "candidate_female"):
        raise ValueError(f"unknown sensory input: {input_name}")
    inputs = network.groups[input_name]
    rate = SENSORY_TOTAL_HZ / len(inputs)
    events = np.zeros((steps, len(inputs)), dtype=np.bool_)
    events[stim_start:stim_end] = (
        np.random.RandomState(seed).random_sample((stim_end - stim_start, len(inputs)))
        < rate
        * DT_MS
        / 1000.0
    )
    return inputs, events


def event_digest(events: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(events).tobytes()).hexdigest()


def run_trial(
    network: Network,
    input_name: str,
    seed: int,
    intervention: Intervention | None = None,
    reversal_mV: float = DEFAULT_REVERSAL_MV,
    signs: np.ndarray | None = None,
    events: np.ndarray | None = None,
) -> TrialResult:
    interv = intervention if intervention is not None else none()
    use_signs = apply_sign_overrides(network) if signs is None else signs
    inputs, generated = sensory_event_schedule(network, input_name, seed)
    use_events = generated if events is None else events
    if events is not None and input_name != "no_input" and use_events.shape[1] != len(inputs):
        raise ValueError("provided events shape does not match input population")
    packed = interv.packed(network.n)
    act_target = interv.target if interv.kind == "activation" else np.zeros(0, np.int32)
    # Class Poisson activation uses a dedicated stream so it does not consume
    # the matched sensory schedule RNG (already generated outside the kernel).
    # Seed the kernel RNG only when class activation is requested.
    if interv.kind == "activation" and interv.activation_hz > 0:
        np.random.seed(seed + 1_000_003)
    counts, vmin, vmax = _simulate(
        network.graph.indptr,
        network.graph.indices,
        network.graph.data,
        use_signs,
        inputs if input_name != "no_input" else np.zeros(0, np.int32),
        use_events,
        packed["output_mask"],
        float(packed["output_gain"]),
        packed["tonic"],
        float(packed["activation_hz"]),
        act_target,
        float(reversal_mV),
        DT_MS,
    )
    if vmin < reversal_mV - 1e-8:
        raise AssertionError(f"voltage {vmin} below reversal {reversal_mV}")
    return TrialResult(
        counts=counts,
        voltage_min_mV=float(vmin),
        voltage_max_mV=float(vmax),
        events=int(use_events.sum()),
        event_sha256=event_digest(use_events),
        input_name=input_name,
        seed=int(seed),
        intervention=interv.name,
        reversal_mV=float(reversal_mV),
    )


def response_metrics(
    result: TrialResult, network: Network, intervention: Intervention | None = None
) -> dict:
    """Summarize P1/pC1, mAL, and global activity in the 250 ms window."""
    p1 = network.groups["P1_readout"]
    broad = network.groups["P1_related"]
    mal = network.groups["mAL"]
    window = result.counts[5:]
    p1_spikes = int(window[:, p1].sum())
    n_bins = window.shape[0]
    window_s = n_bins * BIN_MS / 1000.0
    metrics = {
        "input": result.input_name,
        "seed": result.seed,
        "intervention": result.intervention,
        "reversal_mV": result.reversal_mV,
        "events": result.events,
        "event_sha256": result.event_sha256,
        "P1_spikes": p1_spikes,
        "P1_Hz": float(p1_spikes / len(p1) / window_s),
        "broad_P1_spikes": int(window[:, broad].sum()),
        "mAL_spikes": int(window[:, mal].sum()),
        "network_spikes": int(result.counts.sum()),
        "response_window_spikes": int(window.sum()),
        "active_neurons": int(np.any(result.counts, axis=0).sum()),
        "voltage_min_mV": result.voltage_min_mV,
        "voltage_max_mV": result.voltage_max_mV,
    }
    if intervention is not None and intervention.target.size:
        # Guard: baseline sensory interventions must not target P1 readout.
        overlap = np.intersect1d(intervention.target, p1)
        metrics["intervention_target_size"] = int(intervention.target.size)
        metrics["intervention_overlaps_P1"] = int(overlap.size)
    return metrics


def preference_score(male_p1: float, female_p1: float, eps: float = 1e-9) -> float:
    """Modeled response asymmetry in [-1, 1]; not behavioral preference."""
    denom = male_p1 + female_p1 + eps
    return float((male_p1 - female_p1) / denom)
