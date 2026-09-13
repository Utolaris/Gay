# E2 results: male excitatory relay amplification

Female pathway and sign map intact. Tonic 6 mV frozen in
[E2_PROTOCOL.md](E2_PROTOCOL.md) before outcomes. 10 conditions × 3 seeds ×
2 inputs = 60 runs. Raw: `e2-results.json`.

References: WT male 0,0,0 / female 11,8,8; mAL_silence male 4,5,1 / female 18,18,20.

## Grid (P1 spikes by seed)

| Condition | Male P1 | Female P1 | Mean pref | Male>Female all seeds? | Female unchanged? |
|---|---|---|---:|---|---|
| mAL_silence (ref) | 4, 5, 1 | 18, 18, 20 | −0.70 | no | — |
| vAB3ef_tonic6 | 0, 0, 0 | 7, 8, 7 | −1.00 | no | **no** (F↓) |
| vAB3ef_tonic6 + mAL | 3, 4, 3 | 16, 19, 20 | −0.69 | no | no |
| AN03A008_tonic6 | 0, 0, 0 | 11, 8, 8 | −1.00 | no | yes |
| AN03A008_tonic6 + mAL | **1, 2, 1** | 18, 18, 20 | −0.87 | no | yes |
| AN09B002_tonic6 | 0, 0, 0 | 11, 8, 8 | −1.00 | no | yes |
| AN09B002_tonic6 + mAL | 4, 5, 1 | 18, 18, 20 | −0.70 | no | yes (= mAL alone) |
| all4_tonic6 | 0, 0, 0 | 7, 8, 7 | −1.00 | no | no |
| all4_tonic6 + mAL | 4, 2, 1 | 16, 19, 20 | −0.77 | no | no |

`phenotype_C_conditions`: **[]** (empty).

## Answers

### Can 6 mV tonic on male-specific excitatory relays raise male-cue P1 alone?

**No.** Every relay-only condition leaves male at 0 on all seeds. mAL (and
other brakes) still dominate; subthreshold relay depolarization does not
reach the readout.

### Can it push male above female under mAL silence?

**No.**

- Best male counts remain the mAL-alone reference (4, 5, 1).
- AN03A008 + mAL **reduces** male (1, 2, 1) — amplifying this relay is
  counterproductive in-model (likely recruits extra inhibition / saturation).
- AN09B002 + mAL is identical to mAL alone — no measurable contribution.
- vAB3ef + mAL does not beat mAL alone on male and slightly perturbs female.

### Are these relays male-specific in this encoding?

- **AN03A008 / AN09B002:** yes in the weak sense that female P1 is unchanged —
  but they also fail to move male P1.
- **AN09B017e/f (vAB3):** **not isolated.** Tonic on e/f lowers female P1
  (11,8,8 → 7,8,7) and increases mAL spikes, consistent with the strong
  vAB3→mAL (5437 syn) edge: driving vAB3 co-drives the brake.

## Interpretation (model scope)

1. E2 is a **structured null**: weak male excitatory anatomy cannot be rescued
   by 6 mV tonic drive on the best positive-path intermediates.
2. vAB3 is a poor “male amplifier” in this graph because it is **wired into
   mAL**, not primarily into P1 (direct vAB3→P1 = 8 synapses).
3. Combined with E1: neither “cut another brake” nor “boost the known male
   relay” yields phenotype C with female intact at frozen parameters.
4. Remaining informative directions: **E5** (persistent/arousal state that
   might nonlinearly favor sparse male drive), **E7** (input encoding —
   uniform Poisson may erase labeled lines), **E8** (structural rewiring
   sensitivity for vAB3→P1 / male→P1), and **E6** (broader pC1 readouts).

## Claim discipline

| Statement | Status |
|---|---|
| Tonic 6 mV on AN09B017e/f / AN03A008 / AN09B002 creates male>female | Does not hold |
| vAB3 activation is a clean male-specific lever here | Does not hold (co-drives mAL) |
| AN03A008 amplification helps male under mAL block | Opposite (male decreases) |
| Biological impossibility of male-biased preference | Not claimed — only this grid |

## Reproduce

```sh
uv run python run_e2.py .experiment-data
```
