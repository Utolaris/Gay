# S-block protocol: b/c sign sensitivity + hop-3 local interventions

Frozen before outcomes. Seeds 11–13. LIF constants, graph, female excitatory
override, and matched schedules unchanged from baselines. No direct P1
stimulation. Female pathway is never silenced or sign-flipped.

## Block S1 — AN09B017b/c sign sensitivity (assumption test)

Question: is the male hop-2 “brake” an anatomical fact or a consequence of the
Glu→−1 default?

Evidence (MaleCNS v1.0, local tables): consensus/ground_truth NT for b/c is
glutamate; predicted_nt is acetylcholine (conf ≈ 0.94–0.96); receptorType is
empty; only e/f carry Yu-2010 vAB3 synonyms and are forced +1; g also has vAB3
synonyms but stays −1.

Sign variants (applied only to named types; all other signs identical):

| ID | b | c | g |
|---|---:|---:|---:|
| base | −1 | −1 | −1 |
| bc_exc | +1 | +1 | −1 |
| bc_zero | 0 | 0 | −1 |
| bcg_exc | +1 | +1 | +1 |

Each variant × {WT, mAL_output_silence} × inputs {male, female} × seeds 11–13.

Interpretation rule: this is a **model-assumption sensitivity**, not a
preferred biological sign and not a phenotype search.

## Block S2 — Local female-intact interventions (hop-3 motivated)

Anatomical targets from `analyze_bc_sign_hop3.py` (baseline signs):

- Male-specific inhibitory sink: **AN09B017c** (LgLG6→c ≫ LgLG5→c; large
  negative hop-3 mass as t1).
- Male-biased excitatory relay: **AN03A008** (ACh +1; LgLG6→AN03A008 ≫
  LgLG5; direct →P1 present).
- Male-loaded central inhibitory hub: **IN05B011a**.
- Shared mixed-sign hub: **AN09B017d** (large positive hop-3 product under
  −1 sign; E3 full silence may remove disinhibition with the brake).

Frozen conditions (tonic 6 mV as in E2):

| Condition | Interventions |
|---|---|
| WT | identity |
| mAL | mAL output silence |
| c_mAL | silence AN09B017c + mAL |
| bc_mAL | silence b+c + mAL |
| bcd_mAL | silence b+c+d + mAL (prior best, confirm) |
| a008_mAL | tonic 6 mV on AN03A008 + mAL |
| a008_bc_mAL | AN03A008 tonic + silence b+c + mAL |
| a008_c_mAL | AN03A008 tonic + silence c + mAL |
| in11a_mAL | silence IN05B011a + mAL |
| in11a_bc_mAL | silence IN05B011a + silence b+c + mAL |
| d_g25_mAL | output_gain 0.25 on AN09B017d + mAL |
| a008_in11a_mAL | AN03A008 tonic + silence IN05B011a + mAL |

Grid size: 12 × 2 inputs × 3 seeds = 72 runs. S1: 8 × 2 × 3 = 48 runs.
Total 120 runs.

## Prohibited

- Direct P1 stimulation or tonic drive of the readout.
- Cutting, silencing, or sign-flipping the female sensory pathway.
- Global threshold / weight retuning.
- Seed selection or dropping null / female-biased rows.
- Expanding tonic levels or gain ladders after seeing results.

## Success criteria (pre-stated)

Primary: male P1 rises vs mAL-alone while female stays near the mAL-alone
level (or WT, if the condition has no mAL). Phenotype C only if male > female
on all seeds with female intact.

Secondary: sign-flip of b/c alone reverses hop-2/3 mass ranking (anatomical
sensitivity already computed); does it also reverse P1 drive in LIF?

Null / opposite outcomes are reported.
