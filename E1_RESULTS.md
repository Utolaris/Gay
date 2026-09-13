# E1 results: P1 brake lesion map

**Female pathway intact** (sign map and female excitatory override unchanged).
All 15 conditions × 3 seeds × 2 inputs = 90 runs. Protocol frozen in
[E1_PROTOCOL.md](E1_PROTOCOL.md) before outcomes. Raw rows: `e1-results.json`.

WT reference P1 spikes: male 0,0,0; female 11,8,8.

## Primary contrasts (mean ΔP1 vs WT)

| Lesion | n cells | Δ male | Δ female | ΔM−ΔF | male↑ all seeds? | male-specific? |
|---|---:|---:|---:|---:|---|---|
| **mAL_all** | 75 | **+3.33** | **+9.67** | −6.33 | **yes** (4,5,1) | no |
| SIP112m | 8 | +0.33 | +6.33 | −6.00 | no | no |
| SIP113m | 5 | 0 | +3.33 | −3.33 | no | no |
| AN09B017d | 2 | 0 | −1.67 | +1.67 | no | no |
| AN09B017c | 2 | +0.33 | 0 | +0.33 | no | no |
| mAL_m8 / m5b / m1 / m2b / m5c | 6–16 | 0 | 0 | 0 | no | no |
| SIP100m / VES022 / LH004m / AN09B017b | 2–11 | 0 | 0 | 0 | no | no |

## Answers to the E1 questions

### Is mAL special?

**Yes as a distributed population; no as a unique type-level mechanism.**

- Only **full mAL output silence** unmasks a male-cue P1 response on every
  fixed seed (0 → 4, 5, 1). That reproduces the pre-registered baseline.
- **No single mAL subtype** (m8, m5b, m1, m2b, m5c — the largest by P1
  synapse weight) changes male or female P1 at all. The brake is **redundant
  across subtypes**; lesioning one leaves enough mAL output.

### Are other anatomical P1 inhibitors functional gates here?

**Mostly no, under this 250 ms uniform sensory drive.**

- SIP100m (637 syn into P1), VES022 (414), LH004m (350): **zero** ΔP1 on both
  cues. High anatomical weight ≠ required gate on this timescale/encoding.
- SIP112m and SIP113m **raise female** P1 substantially and barely touch male
  → they look like **female-path / shared** brakes, not male gates.
- AN09B017b/c/d (strong male depth-2 *anatomical* inhibitory products) do
  **not** unmask male P1 when silenced. Anatomical product weight over-predicted
  their short-window causal role (direct P1 synapses of b/c/d are only 7/55/7).

### Is any brake male-specific?

**None.** Criterion was Δmale > 0 and Δmale > Δfemale on all seeds with female
intact. Zero lesions qualify. The best male-raising lesion (mAL_all) raises
female ~3× more.

## Interpretation (model scope only)

1. In this model stack, **mAL is a distributed, redundant brake** whose joint
   removal is necessary for male-cue P1 spikes; subtype and non-mAL lesions
   are insufficient.
2. mAL is **not male-specific**: it also strongly suppresses female-cue P1.
   mAL block → phenotype **B** (bidirectional), not **C** (male-biased).
3. E1 does **not** support “find another brake to flip male>female.” The
   informative next move is **E2**: amplify weak male-specific excitatory
   relays (AN09B017e/f, AN03A008, AN09B002) with female pathway intact.
4. Null SIP/VES/LH lesions are real results, not failures — keep them.

## Claim discipline

| Statement | Status |
|---|---|
| Full mAL silence increases male-cue P1 in this model | Holds (all seeds) |
| A single mAL subtype is the brake | Does not hold |
| SIP100m/VES022/LH004m gate P1 in 250 ms sensory window | Does not hold here |
| Any lesion yields male-specific disinhibition with female intact | Does not hold |
| Biological uniqueness of mAL | Not established by this map |

## Reproduce

```sh
uv run python run_e1.py .experiment-data
```
