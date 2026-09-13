# E1 protocol: P1 brake lesion map (female pathway intact)

Frozen **before** running E1 outcomes. Does not modify baseline sign map,
LIF constants, seeds, or female excitatory override.

## Question

Under matched male/female sensory schedules, which inhibitory classes that
project to the eight-cell P1 readout actually gate male-cue vs female-cue P1
activity when silenced **one class at a time**? Is mAL special among brakes?

## Design

- Graph / signs / overrides: identical to pre-registered baselines
  (vAB3e/f excitatory; candidate_female excitatory intact).
- Reversal: −70 mV.
- Seeds: 11, 12, 13.
- Inputs: candidate_male, candidate_female (no no-input; already silent).
- Intervention: `output_silence` on exactly one class per lesion condition.
- **No** direct P1 stimulation, **no** global threshold change, **no** seed
  selection, **no** female-pathway removal.

## Frozen lesion list

Anatomical freeze rule: inhibitory types with ≥338 synapses into the P1
readout in the prepared graph, plus full mAL as the literature reference, plus
the three strongest male→P1 inhibitory intermediates (AN09B017b/c/d) even
though their direct P1 weights are small (they sit on the male 1-hop hub).

| ID | Class | Rationale |
|---|---|---|
| WT | none | identity control |
| mAL_all | all GABA mAL_m* (75 cells) | literature brake / baseline intervention |
| mAL_m8 | type mAL_m8 | largest mAL subtype → P1 |
| mAL_m5b | type mAL_m5b | 2nd mAL subtype → P1 |
| mAL_m1 | type mAL_m1 | 3rd |
| mAL_m2b | type mAL_m2b | 4th |
| mAL_m5c | type mAL_m5c | 5th |
| SIP100m | type SIP100m | top non-mAL inhibitor → P1 |
| VES022 | type VES022 | non-mAL inhibitor |
| SIP112m | type SIP112m | non-mAL inhibitor |
| LH004m | type LH004m | non-mAL inhibitor |
| SIP113m | type SIP113m | non-mAL inhibitor |
| AN09B017b | type AN09B017b | male-path inhibitory intermediate |
| AN09B017c | type AN09B017c | strongest male depth-2 inhibitory product |
| AN09B017d | type AN09B017d | male-path inhibitory intermediate |

All listed conditions run regardless of outcome.

## Metrics

Per (lesion, seed, input): P1 spikes, broad P1 spikes, mAL spikes, network
spikes, event hash.

Per (lesion, seed): ΔP1 = lesion − WT for male and female; preference score.

Primary contrast: **ΔP1_male vs ΔP1_female** under each single-class silence.
A useful male-specific brake would raise male more than female (or raise male
with little female change). A shared brake would move both similarly.

## Interpretation rules

- Report every cell of the grid.
- Do not call a lesion “male-specific” unless Δmale > Δfemale on all seeds
  with female pathway intact.
- Null results (no change) are informative: that class is not required under
  this encoding.
- Phenotype C (male-biased preference without killing female) is **not** an
  E1 success criterion; E1 is a gate map.
