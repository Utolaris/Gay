# E3–E8 protocol freeze (before outcomes)

Female pathway and baseline sign map stay intact unless an experiment
explicitly states a structural sensitivity (E8) or encoding change (E7).
No direct P1 stimulation. Seeds 11, 12, 13. Reversal −70 mV.
Trial-level process parallelism; single-trial kernel unchanged except where
an experiment introduces a labeled dynamics/encoding/structure change.

## E3 — Selective relief of AN09B017b/c/d inhibition

**Question:** is male-cue P1 sparseness dominated by this local Glu brake
rather than mAL?

| ID | Intervention |
|---|---|
| WT / mAL_ref | references |
| b/c/d_sil | silence each type alone |
| bcd_sil | silence b+c+d together |
| bcd_gain0.25 | output_gain 0.25 on b+c+d |
| bcd_sil_mAL / bcd_gain0.25_mAL | same + full mAL silence |

E1 already showed single-type silence ≈ null; E3 adds combined + gain and the
mAL interaction.

## E4 — Shared-hub E/I routing

Shared male∩female 1-hop hubs: **AN05B035**, **IN05B011a**, **IN05B011b**.

| ID | Intervention |
|---|---|
| hub_sil_* | silence each hub type alone |
| hub_all_sil | silence all three |
| hub_all_gain0.5 | output_gain 0.5 on all three |
| hub_all_gain0.5_mAL | + mAL silence |

**Question:** do shared hubs move male and female differentially?

## E5 — Recurrent / arousal state

New labeled dynamics (not in production `_simulate`): multiplicative synaptic
gain state `S(t)` ∈ [1, S_max]:

- `S ← min(S_max, S + α)` on each network spike (counted once per spiking neuron per step)
- `S ← 1 + (S−1)·exp(−dt/τ_s)` each step
- Frozen: `S_max=2.0`, `α=0.02`, `τ_s=50 ms`
- Applied to all outgoing weights (global arousal), **or** male-relay-gated
  variant: `S` multiplies only edges whose presynaptic type ∈ {AN09B017e,f, AN03A008, AN09B002}

| ID | Variant | + mAL? |
|---|---|---|
| arousal_global | global S | no / yes |
| arousal_malerelay | S only on male-relay outgoing | no / yes |

**Question:** does nonlinear state amplify sparse male drive more than female?

## E6 — Readout robustness

Re-run the key contrast set and score **both** eight-cell pC1_4a/b and the
49-cell broad `P1_related` cohort (and report both).

Conditions: WT, mAL_silence, bcd_sil_mAL, arousal_global_mAL (if E5 runs).

## E7 — Input encoding audit

Matched **total** event count ≈ 1124/1073/1060-class totals (use same
Bernoulli rate on a **subset**).

| ID | Encoding |
|---|---|
| full | all cells of the class (baseline) |
| half_random | first half of sorted indices (deterministic) |
| type_only_a | only LgLG6 (male) / LgLG5 (female) |
| type_only_b | only LgLG7 (male) / LgLG8 (female) |

Each × {WT, mAL_silence}.

**Question:** does uniform all-cell Poisson erase labeled-line structure?

## E8 — Structural rewiring sensitivity (hypothetical)

Sensitivity only — **not** a biological claim. Copy graph; multiply weights of
edges AN09B017e/f → P1_readout by `{10, 100}` (from 8 synapses toward a
literature-implied strong direct path). Everything else unchanged.
Conditions: each factor × {WT, mAL_silence}.

## Success / interpretation rules

- Phenotype C only if male > female on all seeds **with female pathway intact**
  (E8 is explicitly a what-if graph).
- Report all grid cells; nulls are results.
- Do not retune levels after seeing outcomes.
