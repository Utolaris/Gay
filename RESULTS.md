# Orientation baseline results

All values below are measured. Nothing was retuned after looking at outcomes.

## Upstream reproduction

`experiment/followup/bounded_routes.py` was re-run from prepared data in this
worktree. All 18 candidate-male/female P1 spike counts, event counts, and event
SHA-256 hashes matched `experiment/followup/summary.json` and
`figures/figure-data.csv`. The overwritten `bounded-route-results.json` was
restored with `git checkout` so upstream files remain unchanged.

## Orientation baselines (this package)

Protocol: [PROTOCOL.md](PROTOCOL.md). Implementation is independent; LIF
constants and sign overrides match the upstream bounded sensory model.

Fixed grid: interventions `{WT, mAL_output_silence}` × inputs
`{candidate_male, candidate_female, no_input}` × seeds `{11,12,13}` = 18 runs
at reversal −70 mV.

### P1 readout spikes (250 ms window)

| Input | Seed | WT | mAL output silence | Δ |
|---|---:|---:|---:|---:|
| male | 11 | 0 | 4 | +4 |
| male | 12 | 0 | 5 | +5 |
| male | 13 | 0 | 1 | +1 |
| female | 11 | 11 | 18 | +7 |
| female | 12 | 8 | 18 | +10 |
| female | 13 | 8 | 20 | +12 |
| no-input | 11–13 | 0 | 0 | 0 |

Male-cue disinhibition reproduces upstream: WT male response is silent; mAL
output silence yields 1–5 spikes on every fixed seed. Female response remains
stronger in both conditions.

Event schedules are identical within each (input, seed) across interventions.
no-input is silent. Voltage minima stay above −70 mV. No intervention targets
the P1 readout (`intervention_overlaps_P1 = 0` on every silenced run).

### Preference score

`(male − female) / (male + female + 1e-9)` on P1 spike counts.

| Intervention | Seed | Male | Female | Score |
|---|---:|---:|---:|---:|
| WT | 11 | 0 | 11 | −1.000 |
| WT | 12 | 0 | 8 | −1.000 |
| WT | 13 | 0 | 8 | −1.000 |
| mAL silence | 11 | 4 | 18 | −0.636 |
| mAL silence | 12 | 5 | 18 | −0.565 |
| mAL silence | 13 | 1 | 20 | −0.905 |

Mean score: WT ≈ −1.000; mAL silence ≈ −0.702.

**Interpretation:** blocking mAL output enables a previously suppressed
male-cue P1 response and moves the modeled asymmetry toward zero. It does
**not** reverse sign to male-biased. The score is a response index under this
model and input encoding, not mate preference or choice probability.

### Checks

- `pytest`: 11 tiny-network unit tests pass (excitation, inhibition, silence,
  gain bounds, activation, matched schedules, silent no-input, score bounds,
  P1-overlap guard, API validation).
- `compare_upstream.py`: orientation WT/silence P1 counts match upstream
  bounded pairs for all 18 male/female seed combinations at −70 mV.
- Original `experiment/` sources and committed follow-up results are unmodified.

## Reproduce

```sh
# from worktree root, with .experiment-data prepared
uv run --project experiment/orientation pytest -q experiment/orientation/tests
uv run --project experiment/orientation python experiment/orientation/run_baselines.py .experiment-data
uv run --project experiment/orientation python experiment/orientation/compare_upstream.py \
  experiment/orientation/baseline-results.json \
  experiment/followup/summary.json
```
