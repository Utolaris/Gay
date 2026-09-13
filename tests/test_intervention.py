"""Unit tests on tiny networks before interpreting full-network orientation results."""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orientation.intervention import (
    activate_tonic,
    none,
    output_gain,
    output_silence,
)
from orientation.simulate import (
    Network,
    preference_score,
    response_metrics,
    run_trial,
    sensory_event_schedule,
)


def _pair_network(sign_pre: int = 1) -> Network:
    # neuron 0 driven sensory cell, neuron 1 silent target
    graph = csr_matrix(
        (np.array([100.0], dtype=np.float32), (np.array([0]), np.array([1]))),
        shape=(2, 2),
    )
    signs = np.array([sign_pre, 1], dtype=np.int8)
    ids = np.array([100, 200], dtype=np.int64)
    groups = {
        "candidate_male": np.array([0], dtype=np.int32),
        "candidate_female": np.array([], dtype=np.int32),
        "mAL": np.array([0], dtype=np.int32),
        "P1_related": np.array([1], dtype=np.int32),
        "P1_readout": np.array([1], dtype=np.int32),
        "vAB3": np.array([], dtype=np.int32),
    }
    return Network(graph=graph, signs=signs, ids=ids, groups=groups)


def test_excitatory_transmission_drives_target():
    net = _pair_network(sign_pre=1)
    result = run_trial(net, "candidate_male", seed=11, intervention=none())
    assert result.counts[5:, 1].sum() > 0
    assert result.events > 0


def test_inhibitory_sign_does_not_drive_target():
    net = _pair_network(sign_pre=-1)
    result = run_trial(net, "candidate_male", seed=11, intervention=none())
    assert result.counts[5:, 1].sum() == 0


def test_output_silence_blocks_target_but_preserves_presynaptic_spikes():
    net = _pair_network(sign_pre=1)
    intact = run_trial(net, "candidate_male", seed=11, intervention=none())
    silenced = run_trial(
        net, "candidate_male", seed=11, intervention=output_silence([0])
    )
    assert intact.events == silenced.events
    assert intact.event_sha256 == silenced.event_sha256
    # same external schedule => identical sensory-cell spikes
    assert np.array_equal(intact.counts[:, 0], silenced.counts[:, 0])
    assert intact.counts[5:, 1].sum() > 0
    assert silenced.counts[5:, 1].sum() == 0


def test_output_gain_is_bounded_and_monotone():
    net = _pair_network(sign_pre=1)
    g0 = run_trial(net, "candidate_male", seed=11, intervention=output_gain([0], 0.0))
    g1 = run_trial(net, "candidate_male", seed=11, intervention=output_gain([0], 1.0))
    assert g0.counts[5:, 1].sum() == 0
    assert g1.counts[5:, 1].sum() > 0
    # gain > max_gain rejected
    with pytest.raises(ValueError):
        output_gain([0], 1.5, max_gain=1.0)


def test_activation_tonic_does_not_require_p1_targeting():
    # Drive the sensory class itself; target remains network-mediated.
    net = _pair_network(sign_pre=1)
    act = activate_tonic([0], tonic_mV=12.0)
    # zero sensory events: only activation should matter if we pass empty schedule
    empty = np.zeros((3000, 1), dtype=np.bool_)
    result = run_trial(
        net,
        "candidate_male",
        seed=11,
        intervention=act,
        events=empty,
    )
    assert result.events == 0
    assert result.counts[:, 0].sum() > 0
    assert result.counts[5:, 1].sum() >= 0


def test_no_input_is_silent_from_rest():
    net = _pair_network(sign_pre=1)
    result = run_trial(net, "no_input", seed=11, intervention=none())
    assert result.events == 0
    assert result.counts.sum() == 0


def test_matched_schedules_across_interventions():
    net = _pair_network(sign_pre=1)
    _, events = sensory_event_schedule(net, "candidate_male", 12)
    a = run_trial(net, "candidate_male", 12, none(), events=events)
    b = run_trial(net, "candidate_male", 12, output_silence([0]), events=events)
    assert a.event_sha256 == b.event_sha256
    assert a.events == b.events


def test_sensory_schedule_seed_is_deterministic_and_input_specific():
    net = _pair_network(sign_pre=1)
    _, e1 = sensory_event_schedule(net, "candidate_male", 11)
    _, e2 = sensory_event_schedule(net, "candidate_male", 11)
    assert np.array_equal(e1, e2)
    # female population empty in tiny net; use male only and no_input
    inp_n, e_n = sensory_event_schedule(net, "no_input", 11)
    assert e_n.shape[1] == 0 and e_n.sum() == 0


def test_preference_score_bounds():
    assert preference_score(0, 0) == 0.0
    assert preference_score(5, 0) == pytest.approx(1.0)
    assert preference_score(0, 5) == pytest.approx(-1.0)
    score = preference_score(3, 1)
    assert -1.0 < score < 1.0
    assert score > 0


def test_response_metrics_and_p1_overlap_guard():
    net = _pair_network(sign_pre=1)
    interv = output_silence([0])
    result = run_trial(net, "candidate_male", 11, interv)
    metrics = response_metrics(result, net, interv)
    assert metrics["intervention_overlaps_P1"] == 0
    assert "P1_spikes" in metrics and "network_spikes" in metrics
    bad = output_silence([1])  # targets P1 in this tiny net
    result_bad = run_trial(net, "candidate_male", 11, bad)
    metrics_bad = response_metrics(result_bad, net, bad)
    assert metrics_bad["intervention_overlaps_P1"] == 1


def test_intervention_validation():
    with pytest.raises(ValueError):
        output_silence([])
    with pytest.raises(ValueError):
        output_gain([0], -0.1)
    with pytest.raises(ValueError):
        activate_tonic([0], 0.0)
