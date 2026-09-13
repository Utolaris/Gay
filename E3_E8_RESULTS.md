# E3–E8 results

Protocol frozen in [E3_E8_PROTOCOL.md](E3_E8_PROTOCOL.md). 186 jobs, seeds
11–13, female pathway intact except E8 structural what-if. Raw:
`e3-e8-results.json`. Broad 49-cell readout is in every row (E6).

References: WT male 0,0,0 / female 11,8,8; mAL_silence male 4,5,1 / female 18,18,20.

**No condition achieves phenotype C** (male > female on all seeds with female
intact). Below: what each axis actually did.

## E3 — AN09B017b/c/d selective relief

| Condition | Male P1 | Female P1 | pref |
|---|---|---|---:|
| b/c/d silence alone | ≈0–1 | ≈6–11 | ≈−1.0 |
| **bcd_sil + mAL** | **9, 9, 5** | 18, 18, 20 | **−0.42** |
| bcd_gain0.25 + mAL | 7, 4, 5 | 16, 18, 20 | −0.54 |
| mAL alone | 4, 5, 1 | 18, 18, 20 | −0.70 |

**Finding:** the local Glu brake **does** gate male once mAL is off. Combined
bcd silence + mAL is the **best male drive so far with female intact**
(9,9,5). Still female-biased. Single-type silence remains weak (matches E1).

## E4 — Shared hubs (AN05B035, IN05B011a/b)

All hub lesions leave male at 0 (or worse under mAL). hubs_gain0.5 + mAL →
male 4,0,0. Shared hubs are not a male-specific lever in this encoding.

## E5 — Arousal state (frozen S_max=2, α=0.02, τ=50 ms)

| Condition | Male | Female | pref |
|---|---|---|---:|
| arousal_global | 2,2,3 | 39,44,34 | −0.88 |
| **arousal_global + mAL** | **29,27,38** | **68,57,56** | **−0.32** |
| arousal on male-relay only | 0 | =WT | −1.0 |
| male-relay + mAL | = mAL alone | = mAL | −0.70 |

**Finding:** global arousal strongly amplifies **both** sexes once mAL is off.
Male reaches 27–38 spikes (highest absolute male so far) but female remains
~2× higher. Relay-gated arousal is inert (relays barely spike).

## E6 — Broad readout

Every summary row carries `male_broad` / `female_broad`. Qualitative ranking
of conditions matches the 8-cell P1 readout; broad cohort does not reverse
male vs female.

## E7 — Input encoding

| Encoding | Male under mAL | Female under mAL |
|---|---|---|
| full | 4,5,1 | 18,18,20 |
| half cells | 2,0,0 | 21,22,25 |
| type_only_a (LgLG6 / LgLG5) | 1,0,1 | 22,17,19 |
| **type_only_b (LgLG7 / LgLG8)** | **0,0,0** | **0,0,0** |

**Finding:** uniform Poisson is not erasing a strong male labeled line — male
is weak under every encoding. LgLG8-only female drive is **silent** even with
mAL block; LgLG5 carries the female drive. Encoding changes do not create
male bias.

## E8 — vAB3e/f→P1 rewiring sensitivity (hypothetical graph)

| Weights × | Condition | Male | Female | pref |
|---:|---|---|---|---:|
| 10 | WT / +mAL | 0 / 6–8 | ~9 / ~20 | −1.0 / −0.50 |
| **100** | WT | **31,33,30** | **51,52,55** | **−0.25** |
| **100** | +mAL | **41,42,40** | **58,56,61** | **−0.17** |

**Finding:** restoring a literature-scale direct vAB3→P1 path (from 8
synapses ×100) produces large male responses and moves preference score
closest to zero **without deleting the female pathway**. Female still wins.
This is a **mapping-discrepancy sensitivity**, not a biological result.

## Cross-experiment ranking (female intact, mean preference)

| Rank | Condition | Male | Female | pref |
|---:|---|---|---|---:|
| 1 | e8×100 + mAL | 40–42 | 56–61 | −0.17 |
| 2 | e8×100 WT | 30–33 | 51–55 | −0.25 |
| 3 | arousal_global + mAL | 27–38 | 56–68 | −0.32 |
| 4 | bcd_sil + mAL | 5–9 | 18–20 | −0.42 |

## Interpretation (model scope)

1. **E3 is a real mechanistic hit:** AN09B017b/c/d jointly suppress male-cue
   P1 when mAL is also silenced — complementary brakes, not interchangeable
   single types.
2. **E5** shows sparse male drive *can* be amplified, but global state
   amplifies female more (shared hubs).
3. **E7** argues against “uniform encoding hides the male line.”
4. **E8** identifies the vAB3→P1 **anatomical deficit** as the structural
   bottleneck closest to a male-biased (or at least balanced) readout —
   without touching female signs.
5. Still **no phenotype C** under any female-intact condition.

## Claim discipline

| Statement | Status |
|---|---|
| b/c/d + mAL raises male more than mAL alone | Holds (9,9,5 vs 4,5,1) |
| Male-specific brake exists as a single type | Does not hold |
| Arousal creates male bias | Does not hold (both sexes up) |
| Uniform encoding is the sole cause of male silence | Does not hold |
| vAB3→P1 deficit is the main structural gap for male drive | Supported as **sensitivity**, not biology |

## Reproduce

```sh
uv run python run_e3_e8.py .experiment-data --workers 6
```
