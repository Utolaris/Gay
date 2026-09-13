# Exploratory male-bias parameter search (POST-HOC)

**This is not the pre-registered baseline.** `PROTOCOL.md` and
`baseline-results.json` are unchanged. Every condition below was chosen after
seeing that WT/mAL-silence stays female-biased. Do not present these numbers as
confirmatory evidence of male preference.

Full machine-readable output: `exploratory-male-bias.json`.

Fixed across all exploratory runs: seeds 11/12/13, reversal −70 mV, full graph,
no direct P1 stimulation, matched male/female schedules within each seed.

## Question

Which interventions make male-cue P1 spike count exceed female-cue P1 count on
every fixed seed?

## What did **not** flip bias

| Lever | Range tried | Result |
|---|---|---|
| Male sensory total rate | 5550 → 11100 / 22200 / 44400 Hz | Male P1 stays ~0–6; female still 18–20 under mAL silence |
| Male-class tonic drive | 5 / 10 / 14 mV | No male WT response; under mAL silence male even drops in some seeds |
| Female sensory rate | 5550 → 2775 / 1387.5 Hz | Female still dominates (≈13–19 vs male ≈1–5) |
| Female output gain | 1.0 → 0.5 | Still female-biased |
| Female output gain | 1.0 → 0.25 + mAL silence | Mean still female (3.3 vs 4.7); only 1/3 seeds male>female |

Strengthening the male cue alone does not create a male-biased P1 readout in
this model. The male→P1 path under mAL block remains sparse (1–5 spikes).

## What **did** flip bias on all three seeds

| Condition | Male P1 (seeds) | Female P1 | Preference | Mechanism |
|---|---|---|---|---|
| `no_fem_exc_mAL` | 4, 5, 1 | 0, 0, 0 | ≈ +1.0 | Drop candidate_female excitatory sign override (revert to conservative NT signs) + silence mAL output |
| `fem_gain0.0_mAL` | 4, 5, 1 | 0, 0, 0 | ≈ +1.0 | Force female sensory class output gain to 0 + silence mAL |
| `combo_noexc_male_x2_mAL` | 3, 3, 1 | 0, 0, 0 | ≈ +1.0 | Same as first, with male rate ×2 (does not help male absolute count) |

Reference (unchanged) `ref_mAL_silence` for contrast: male 4, 5, 1; female 18,
18, 20; preference ≈ −0.70.

## Interpretation

1. **The female-biased baseline is driven by the candidate_female excitatory
   override.** LgLG5/8 are glutamate/unclear; forcing them to +1 is what lets
   female cues drive P1. Reverting that override silences the female pathway
   entirely under conservative signs.
2. **mAL output silence is what unmasks any male-cue P1 response.** Without it,
   male is 0 even when female is also 0 (`no_fem_exc_WT`).
3. **Male bias in this model is obtained by removing/weakening the female
   pathway, not by amplifying the male pathway.** That is a statement about
   this graph encoding and sign map, not about biological mate preference.
4. These are post-hoc assumption changes. They cannot replace the pre-registered
   result: under the original matched sign map, mAL block increases male
   response but leaves a female-biased readout.

## Reproduce

```sh
uv run python explore_male_bias.py
```

Requires `.experiment-data/` prepared as in the main README.
