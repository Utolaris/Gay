# REWARD learning results

Protocol frozen in [REWARD_PROTOCOL.md](REWARD_PROTOCOL.md). 18 jobs
(6 conditions × seeds 11–13), 20 training epochs + 20 reversal epochs each.
Raw: `reward-results.json` (`n_ok=18`, `n_error=0`). η=0.05, w_floor=0.25,
τ_e=50 ms, PAM 20 Hz, seeds fixed. No mAL silence, no b/c sign change, no P1 drive.

**No condition produces phenotype C.** Preference stays female-biased
(`post_pref = −1.0` on all seeds). Below: what failed and why.

---

## 0. Anatomical gate (measured before learning)

| Link | Synapses |
|---|---:|
| LgLG → KC | **0** |
| KC → MBON | 61,210 edges / 463,640 syn |
| MBON → P1 | **0** |
| MBON → P1 hop-2 | present (OA-VPM3, SIP106m, …) |

Uniform LgLG Poisson (300 ms) evokes **KC spikes = 0** on male and female
for seeds 11 and 13. Seed 12 male evokes 863 KC spikes, but only in
**bins ~260–300 ms** — after the frozen DA/sensory gate (50–250 ms).

So the MB is anatomically disconnected from the pheromone GRNs at hop-1,
and functionally almost silent or late under this encoding.

---

## 1. Condition summary (post-test P1, seeds 11/12/13)

| Condition | post P1 M | post P1 F | Δw (mean ratio) | KC M pre | MBON M pre→post |
|---|---|---|---:|---|---|
| `male_rew` | 0,0,0 | 11,8,8 | **1.000** (no move) | 0,863,0 | 3/69/3 → same |
| `no_rew` | 0,0,0 | 11,8,8 | 1.000 | 0,863,0 | same |
| `rand_rew` | 0,0,0 | 11,8,8 | 1.000 | 0,863,0 | same |
| `fem_rew` | 0,0,0 | 11,8,8 | 1.000 | 0,863,0 | same |
| `plast_off` | 0,0,0 | 11,8,8 | 1.000 | 0,863,0 | same |
| **`mb_access`** | **0,0,0** | 11,8,8 | **0.92–0.96** (moved) | 0,863,0 | ~same at test |

Reversal training (female+PAM, 20 epochs) does not change any post2 P1.

---

## 2. Why `male_rew` did not learn

Weights never left init (`frac_depressed = 0`) despite seed-12 KC spikes.

Mechanism: those spikes fall **outside the DA window**. Eligibility is set
only when a KC spikes; DA depression applies only on steps with `da_flag`
(50–250 ms). Late KC activity (≈260 ms+) never meets the gate.

This is **temporal misalignment**, not a broken three-factor rule —
`mb_access` uses the same rule and does depress weights.

---

## 3. What `mb_access` shows (positive control)

Pre-registered stand-in for the missing LgLG→KC edge: Poisson 20 Hz on 200
KCs during male cue + PAM reward.

| During training (male trial) | Result |
|---|---|
| KC spikes | ~950–1020 (forced) |
| MBON spikes | 12–27 (up from ~3 at rest) |
| PAM spikes | ~1600 |
| KC→MBON weights | mean ratio → 0.92–0.96; ~5–11% edges at floor |
| **P1 male** | **0 every epoch, every seed** |
| P1 female (unpaired) | 8–11 (unchanged) |

Post-test (DA off, no KC drive): weights remain depressed, but male P1 is
still 0 and MBON returns to baseline. **Association is written into KC→MBON
but never read out**, because the cue still does not drive KCs at test.

Even *during* training, when KCs and PAM are both active, P1 stays 0:
MBON activity is insufficient (or sign-routed) to recruit the courtship
readout under frozen b/c signs and intact mAL.

---

## 4. Answers to the three core questions

### Q1 — Does male cue enter the plastic circuit?

**Mostly no.** Direct LgLG→KC = 0. Evoked KC activity is 0 (seeds 11,13)
or late and outside the reward gate (seed 12). The MB is not a reliable
readout of the pheromone GRN channel in this model stack.

### Q2 — Can DA conditioning alone form persistent male bias?

**No**, under frozen signs, intact mAL, no P1 drive. `male_rew` does not
move weights and does not move P1. Controls (`no_rew`, `rand_rew`,
`fem_rew`, `plast_off`) are identical nulls on P1.

### Q3 — Learning failure vs innate brake?

**Both, in sequence:**

1. **Learning failure at the access stage** — cue does not co-activate KCs
   with DA, so the three-factor rule has nothing to bind.
2. **Expression failure even when association is forced** — `mb_access`
   depresses KC→MBON and raises MBON spikes during training, yet P1 remains
   0. Residual drive is blocked before or at the courtship readout
   (mAL / hop-2 inhibitory routing / weak MBON→P1), **not** by b/c sign
   (unchanged here) and not by plasticity off.

Distinguishing labels used in the protocol:

| Label | Holds? |
|---|---|
| Reward association learned (weights/MBON change) | Only `mb_access` |
| Courtship output expressed (P1 male bias) | **Never** |
| Innate brake after association | **Yes** (mb_access) |

---

## 5. Mechanism conclusions

1. **Mushroom-body reward learning is not a viable path to male-biased P1
   in this connectome + encoding** without adding a sensory→KC access edge
   (or a longer/stronger olfactory channel that this package does not model).
2. **Dopamine as a plasticity gate is necessary but not sufficient.** It
   only affects active KC–MBON synapses; silent or late KCs yield no
   learning (all `male_rew` controls).
3. **Even a successful KC→MBON depression does not release male P1** while
   mAL is intact and b/c remain inhibitory. The innate brake story from
   S-block stands: expression is gated downstream of MB.
4. **Do not “fix” this by flipping b/c or silencing mAL inside the learning
   block** — those are separate interventions already measured. This block’s
   job was to test reward-only plasticity; the answer is negative.

---

## 6. Claim discipline

| Statement | Status |
|---|---|
| LgLG does not innervate KC in MaleCNS v1.0 | Anatomical fact |
| Uniform male cue does not drive KC in the DA window | Model result (3 seeds; seed12 late only) |
| DA-gated KC→MBON depression works when KCs are co-driven | Model result (`mb_access`) |
| That depression raises male P1 | Does not hold |
| Reward learning produces mate preference | Not claimed |
| MB is unnecessary for courtship in vivo | Not claimed |

---

## 7. Reproduce

```sh
uv run python run_reward.py .experiment-data --workers 6 --out reward-results.json
```
