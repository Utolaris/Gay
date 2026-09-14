# REWARD protocol: male cue + DA-gated local plasticity

Frozen **before** learning outcomes. Topology, LIF constants, b/c signs,
female excitatory override, and seeds unchanged from package baselines.
No mAL silence. No P1 stimulation. No global plasticity.

## Anatomical premises (MaleCNS v1.0)

| Fact | Measurement |
|---|---|
| KC population | 4,053 cells (KCg-m/ab/a′b′ classes) |
| MBON | 97 cells |
| PAM | 316 cells |
| LgLG → KC direct | **0 synapses** |
| KC → MBON | 61,210 edges, 463,640 synapses |
| MBON → P1 direct | **0 synapses** |
| MBON → P1 hop-2 | via OA-VPM3, AOTU100m, SIP106m, SMP163, … |
| LgLG → KC within 4 hops | all KCs reachable (hop 3–4), **not** hop-1 |
| Cue-evoked KC spikes (WT 300 ms) | **male 0, female 0** |

Implication: the MB is present but **not driven** by the uniform LgLG encoding
on this timescale. Plasticity cannot form a cue–reward association unless KCs
are active. This is pre-registered as a likely primary null.

## Plasticity site

Only **KC → MBON** anatomical edges are plastic. Three-factor rule
(mushroom-body-canonical dopamine depression of active KC–MBON synapses):

```
elig_i ← 1 on KC spike i
elig_i ← elig_i · exp(−dt / τ_e) each step
if DA(t) > 0:
  w_ij ← clip(w_ij − η · DA(t) · elig_i, w_floor · w0_ij, w0_ij)
```

- `τ_e = 50 ms`, `η = 0.05`, `DA = 1.0` when reward channel is on else 0
- `w_floor = 0.25` (synapses may depress to 25% of anatomical init, never
  reverse sign, never exceed init)
- Reward = class Poisson **20 Hz on PAM** during the sensory window
  (50–250 ms) of **paired trials only**. PAM is never a P1 target.
- No plasticity outside KC→MBON. No weight changes in post-test.

## Phases

1. **Pre-test** (plasticity off, DA off): male and female cue, seeds 11–13.
2. **Training**: `n_epochs = 20` paired trials (male cue + PAM reward),
   then one unpaired female trial (no reward) per epoch. Plasticity on.
   Same seed stream per (condition, seed); schedules matched across
   conditions via identical `RandomState(seed)` sensory generation.
3. **Post-test** (plasticity off, DA off): male and female cue, same seeds.
4. **Reversal** (only if a male bias metric moves; still pre-registered to
   run regardless of direction): 20 epochs female+PAM, then post-test-2.

## Conditions (all × seeds 11,12,13)

| ID | Training pairing | Plasticity |
|---|---|---|
| `male_rew` | male + PAM; female no-reward | on |
| `no_rew` | male no-reward; female no-reward | on |
| `rand_rew` | 50% of male trials get PAM (coin by seed) | on |
| `fem_rew` | female + PAM; male no-reward | on |
| `plast_off` | male + PAM | **off** (DA present, no Δw) |
| `mb_access` | male + PAM **and** Poisson 20 Hz on 200 KCs | on |

`mb_access` is a **pre-registered positive control** for the missing
LgLG→KC edge (stand-in for unmodeled olfactory projection). It does not
modify connectome signs or P1. If `male_rew` is null and `mb_access` is not,
the bottleneck is sensory access to MB, not the DA rule. If `mb_access`
learns MBON changes but P1 stays flat, expression is blocked downstream
(innate brake), not learning.

## Measurements

Per trial: P1 spikes, broad pC1 cohort, KC spikes, MBON spikes, PAM spikes,
network spikes, schedule hash.

Per training epoch: mean |Δw|, fraction of depressed plastic edges, total
plastic weight mass.

Per post-test: ΔP1 male vs female vs pre-test; preference index
`(M−F)/(M+F+ε)`; MBON activity change; residual |Δw| (should be frozen).

Phenotype labels post hoc: learned association (MBON/weight change),
expressed bias (P1 male>female all seeds with female path intact),
blocked expression (weights/MBON change, P1 unchanged).

## Prohibited

- Changing b/c signs, silencing mAL, stimulating P1, global threshold retune.
- Plasticity outside KC→MBON.
- Reward delivery as direct depolarizing current onto P1 or courtship-only
  classes (PAM only).
- Post-hoc η / epoch / floor changes after seeing outcomes.
- Dropping null or opposite-direction runs.

## Core questions

1. Does male cue enter the plastic circuit? (KC activity under cue)
2. Can DA conditioning alone form persistent male bias without sign/mAL edits?
3. If association appears but P1 does not, is innate inhibition blocking
   expression?
