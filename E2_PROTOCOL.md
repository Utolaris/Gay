# E2 protocol: amplify male-specific excitatory relays

Frozen **before** E2 outcomes. Female pathway, sign map, LIF constants, and
seeds unchanged from baselines. No direct P1 stimulation.

## Rationale (from pathway decomposition + E1)

Male→P1 depth-2 excitatory anatomical products are weak compared with
inhibitory AN09B017b/c/d. E1 showed no single brake lesion yields
male-specific disinhibition with female intact. E2 asks whether **amplifying
weak male excitatory relays** can raise male-cue P1 without touching the
female pathway.

## Frozen targets

Types with positive male depth-2 path_sign into P1 (from `analyze_pathways.py`):

| ID | Types | Role |
|---|---|---|
| vAB3ef | AN09B017e + AN09B017f | literature vAB3 pair; forced excitatory |
| AN03A008 | AN03A008 | ACh excitatory intermediate |
| AN09B002 | AN09B002 | ACh excitatory intermediate |
| all4 | e+f+AN03A008+AN09B002 | combined relay set |

## Frozen interventions

- Tonic steady-state depolarizing drive **6 mV** on the relay class only
  (between the 7 mV rest–threshold gap and the exploratory 10 mV; chosen
  analytically, not fitted).
- Two brake states: `mAL_intact`, `mAL_silence`.
- References: `WT`, `mAL_silence` alone (already measured; re-run for pairing).

Grid: 4 relay sets × 2 brake states + 2 references = **10 conditions**.
Seeds 11, 12, 13. Inputs: candidate_male, candidate_female.
Reversal −70 mV. Matched schedules within each (input, seed).

## Success criteria (pre-stated)

- **Primary:** male P1 under (relay + mAL_silence) vs female P1 under the same
  condition, all seeds. Phenotype C only if male > female **with female
  pathway intact** on all seeds.
- **Secondary:** Δmale vs WT / vs mAL_silence alone; whether relay activation
  changes female (it should not, if truly male-specific in this encoding).
- Amplification of female or global runaway is a negative/information result.

## Non-goals

- No preference-score tuning after the fact.
- No changing tonic level if the first grid is female-biased.
- No claiming biological necessity of these types.
