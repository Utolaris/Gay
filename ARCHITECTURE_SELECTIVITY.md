# Architecture synthesis: is sex selectivity = mAL + AN09B017c + FLA?

Question: does real MaleCNS form courtship-circuit sex selectivity by
combining (i) stable male-specific mAL inhibition, (ii) AN09B017c
target-specific action, and (iii) female-biased FLA excitation?

Evidence from this package only (delivered-current LIF on MaleCNS v1.0).
Not a behavioral claim.

---

## 1. Verdict

**Partially yes, with two hard caveats.**

| Component | Role in this model | Sign-robust? | Anatomical? |
|---|---|---|---|
| **mAL series** (esp. **mAL_m8**, mAL_m2b) | Male-only **I** onto P1 (I_f=0 for m8) | **Yes** (A/B/C) | GABA mAL; strong →P1 |
| **AN09B017c** | Male-biased **loading** (LgLG6→c ≫ LgLG5→c) | Loading yes; **effect sign-dependent** | Yes (connectome) |
| **FLA001m / FLA003m** | **Female-biased E** onto P1 (ΔE ≈ −300 to −400) | **Yes** (A/B/C) | Strong →P1 |

Together they form a **last-hop E/I population code** that, under the
package’s Glu→−1 map, yields female-biased P1. They do **not** form a
single labeled-line “sex detector,” and the mechanism is **not robust** to
the b/c transmitter/sign uncertainty.

---

## 2. What each piece actually does

### 2.1 mAL series — sign-robust male brake

Delivered inhibitory flux to P1 (male / female, mean weight-units):

| Type | Map A | Map B | Map C |
|---|---|---|---|
| **mAL_m8** | 582 / 0 | 533 / 0 | 579 / 0 |
| **mAL_m2b** | 246 / 55 | 304 / 55 | 209 / 55 |

- mAL_m8 is the **#1 male-specific inhibitory input** on every sign map.
- Cutting mAL_m8 alone: male rises on B and C; **not enough on A**.
- Full mAL class silence raises **both** sexes (prior E1/S2) — mAL is a
  large shared brake with **subtype-level male bias at m8**.

**This part of the hypothesis holds** as a model/architecture statement.

### 2.2 AN09B017c — target-specific wiring, not a robust gate

| Fact | Status |
|---|---|
| LgLG6→c ≫ LgLG5→c (~170×/cell) | Anatomical |
| c→P1 exists (55 syn) | Anatomical |
| c delivers I to P1 only under male cue | **Only if c=−1** |
| Same edge as E if c=+1 | Map B |
| Silent if c=0 | Map C |

So c is a **male-biased hub edge** (target-specific routing). Calling it a
“male-specific inhibitory gate” **depends entirely on Glu→−1**. Under the
ACh-predicted / vAB3-class map it becomes male-biased **excitation**, and
WT is already male>female.

**Hypothesis clause (ii) is true for connectivity, false as a sign-robust
inhibitory mechanism.**

### 2.3 FLA — female-biased excitation

| Type | E_male / E_female (A) | B | C |
|---|---|---|---|
| FLA001m | 1020 / **1330** | 1092 / **1330** | 1019 / **1330** |
| FLA003m | 602 / **1009** | 568 / **1009** | 608 / **1009** |

Female cue delivers **more excitation** to P1 via FLA on all maps. This is
a **parallel female advantage**, independent of b/c signs, and is **not**
removed by male-specific I cut-sets (explains why female stays at WT when
m8+c are silenced).

**Hypothesis clause (iii) holds** in this model.

---

## 3. How they combine (architecture)

```text
Male cue ──► LgLG6 ──► AN09B017 hub ──► mAL_m8 / m2b / VES022 / PVLP048
                      │                      │
                      │ (c load, sign?)      ▼
                      └──────────────► P1  ◄── FLA001m/003m (female-biased E)
                                         ▲
                                    E/I veto
```

Sex selectivity at the readout is:

1. **Male-specific I list** (mAL_m8 first) that male cue injects into P1
   and female cue does not;
2. **Female-biased E list** (FLA*) that female cue injects more strongly;
3. Optional **c contribution** whose **sign** decides whether the male hub
   adds I (A) or E (B) or nothing (C).

Under Map A this stack is female-biased. Under Map B the same anatomy is
already male-biased without any cut. **The combination does not uniquely
determine sex selectivity — the sign map does.**

---

## 4. What “共同形成性别选择” would require

| Criterion | Met? |
|---|---|
| Three components all present in connectome | **Yes** |
| Each has a causal last-hop effect on P1 in the model | **Yes** (with c sign caveat) |
| Together they explain female>male baseline | **Yes under Map A** |
| Mechanism stable under transmitter uncertainty | **No** (Map B flips) |
| Small cut-set releases male without touching female | **Partly** (m8+c on A; m8 alone on B/C) |
| In vivo / behavioral sex selectivity | **Not tested here** |

---

## 5. Final answer

**In MaleCNS + this LIF stack: yes as a three-part last-hop code; no as a
sign-robust unique mechanism.**

- **Stable male-specific mAL inhibition** — supported (mAL_m8 especially).
- **AN09B017c target-specificity** — supported as **wiring**; its
  **inhibitory gate** identity is **not** robust.
- **FLA female-biased excitation** — supported across sign maps.

These jointly **can** produce sex-biased P1, but the polarity of the male
hub (c) is the fragile hinge. Real MaleCNS therefore supports:

> **Distributed last-hop E/I asymmetry** (male-biased I + female-biased E),
> with **one high-leverage male I subtype (mAL_m8)** and a **sign-unstable
> AN09B017c contribution** —

not a single hard-wired sex switch, and not something this package can
upgrade to a claim about live-fly mate choice.

---

## 6. Claim discipline

| Allowed | Not allowed |
|---|---|
| mAL_m8 is a sign-robust male-specific P1 I source in-model | mAL_m8 is the unique biological sex gate |
| c loading is male-biased; c effect sign is uncertain | c is proven inhibitory in vivo |
| FLA delivers more E to P1 under female cue | FLA causes female preference in behavior |
| Three-part last-hop code explains model sex bias under Map A | This is the in vivo mechanism of sexual preference |

## 7. Key artifacts

| Block | File |
|---|---|
| Delivered-current autopsy | `P1_CURRENT_REPORT.md` |
| Sign robustness | `SIGN_ROBUST_REPORT.md` |
| Pathway death point | `MALE_SIGNAL_AUTOPSY.md` |
| b/c sign adjudication | `BC_SIGN_ADJUDICATION.md` |
