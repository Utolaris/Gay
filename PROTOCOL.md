# Orientation intervention protocol (fixed before outcomes)

Question: under the same bounded full-connectome LIF assumptions as the upstream follow-up, do class-level output interventions on mAL (or other non-readout classes) change P1/pC1 responses and a matched male-vs-female response score?

This directory is independent of `experiment/` and `experiment/followup/` implementation files. It does not modify the connectome graph, LIF constants, or committed upstream result artifacts.

## Fixed model (copied assumptions, not retuned)

- Graph: prepared MaleCNS v1.0 classified-neuron graph (166,606 neurons; 25,574,615 edges).
- Dynamics: bounded inhibitory conductance with finite reversal, matching `experiment/followup/bounded_routes.py` equations.
- Constants: dt 0.1 ms; delay 1.8 ms; refractory 2.2 ms; rest/reset −52 mV; threshold −45 mV; membrane τ 20 ms; synaptic decay 5 ms; anatomical weight scale 0.275; external impulse 68.75 mV; stimulated cells zero refractory; refractory arrivals discarded.
- Primary inhibitory reversal: −70 mV. Optional sensitivity −80 mV is reported only if the baseline grid is expanded; not used to select outcomes.
- Fast signs: ACh +1, GABA −1, glutamate −1, other/unclear 0. Same upstream overrides present in every sensory run: vAB3e/f bodies 11998, 13341, 13693, 512498 set excitatory; candidate_female set excitatory. These are model assumptions, not measured receptor signs.
- Duration 300 ms; sensory window 50–250 ms; readout window bins 5–30 (10 ms bins).

## Populations

- Male cue: `candidate_male` (LgLG6/7 ProLN), 37 cells.
- Female cue: `candidate_female` (LgLG5/8 ProLN), 27 cells.
- No input: empty population, identical timing scaffolding and zero events.
- mAL class: 75 GABA `mAL_m*` cells.
- P1/pC1 readout: eight literature-mapped body IDs 12442, 16719, 17867, 20117, 20803, 23968, 519518, 522419 (pC1_4a/b). Not directly stimulated in sensory baselines.
- Secondary readout retained in metrics only: original 49-cell `P1_related` cohort.

## Intervention API

A neuron-class intervention is one of:

1. `output_silence` — suppress outgoing synaptic transmission from the class (cells may still spike).
2. `activation` — tonic steady-state depolarizing voltage contribution and/or class-specific Poisson external drive. Used for mechanism tests and future interventions, **not** applied to the P1 readout in baseline sensory runs.
3. `output_gain` — multiply outgoing weights of the class by a bounded gain in `[0, max_gain]` with `max_gain = 1.0` by default. `gain=0` is equivalent to output silence; `gain=1` is identity.

Global firing threshold, resting potential, and anatomical weights are not intervention targets.

## Matched conditions

For each seed and each of male / female / no-input:

- External event schedule is generated **once** per (input, seed) with `numpy.random.RandomState(seed)` and reused for every intervention in that input condition.
- Graph, signs, delays, refractory handling, and readout indices are identical across interventions.
- Intervention pairs are therefore only differing by the named class operation.

## Fixed baseline grid

- Interventions: `WT` (identity), `mAL_output_silence`.
- Inputs: male, female, no-input.
- Seeds: 11, 12, 13 (same as upstream bounded follow-up).
- Reversal: −70 mV.
- Total: 2 × 3 × 3 = 18 runs. All runs are executed regardless of outcome.

## Measurements

Per run:

- P1 spikes and Hz in the 250 ms response window.
- Broad P1-related spikes (49-cell cohort).
- mAL spikes.
- Global network spikes and active-neuron count.
- External event count and schedule SHA-256.
- Voltage minimum (assert ≥ reversal − 1e-8).

Per matched (intervention, seed) across male vs female:

- Preference score = `(M − F) / (M + F + 1e-9)` using P1 spike counts in the response window.
- If both M and F are zero, score is 0.

The score is a **modeled response asymmetry index**, not behavioral preference, choice probability, or a fit to courtship.

## Prohibited shortcuts

Not used to obtain a target outcome:

- Direct stimulation or tonic drive of the P1/pC1 readout in sensory baselines.
- Seed selection after looking at results; seeds 11–13 are fixed in advance.
- Changing global threshold, resting potential, anatomical weights, or sign map after seeing outcomes.
- Dropping null or opposite-direction runs.

## Tests before interpreting full-network results

Tiny-network unit tests must show:

- Excitatory transmission drives a silent target.
- Inhibitory sign does not.
- `output_silence` blocks transmission without changing presynaptic spikes under the same external schedule.
- `output_gain` scales outgoing effect monotonically between 0 and 1.
- `activation` tonic drive depolarizes the class without stimulating other classes directly.
- Paired WT vs intervention share identical external schedules.
- No-input is silent from rest.

Full-network checks:

- Matched event hashes within each input/seed across interventions.
- No P1 membership in any baseline activation target.
- Voltage bounds hold.

## Interpretation rule

Report all 18 runs. An increase in male-cue P1 spikes under mAL output silence is conditional model evidence of disinhibition, not male preference. A nonzero preference score is a response ratio artifact of this model and input encoding unless independently validated.

If WT male response is already nonzero, or mAL silence does not increase male response, that outcome is retained and reported without retuning.
