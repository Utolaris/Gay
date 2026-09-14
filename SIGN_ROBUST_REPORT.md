# SIGN_ROBUST_CUTSET results

Protocol: [SIGN_ROBUST_PROTOCOL.md](SIGN_ROBUST_PROTOCOL.md).
Raw: `sign-robust-results.json`. Seeds 11–13. Female path intact.
Each map recomputes delivered-current rankings from scratch.

---

## 0. WT baseline per map (no cut)

| Map | b/c | Male P1 | Female P1 | Notes |
|---|---|---|---|---|
| **A_base** | −1 | **0,0,0** | 11,8,8 | package baseline |
| **B_exc** | +1 | **14,11,11** | 11,8,8 | **male > female already** |
| **C_zero** | 0 | **0,0,1** | 11,8,8 | uncertainty control |

Map B alone produces phenotype C without any ablation (sign sensitivity).

---

## 1. Is mAL_m8 a male-specific brake on all maps?

**Yes.** Delivered I to P1:

| Map | I_male | I_female | score rank |
|---|---:|---:|---|
| A | 582 | 0 | #1 |
| B | 533 | 0 | #1 |
| C | 579 | 0 | #1 |

mAL_m8 is **sign-robust** as a male-specific last-hop inhibitory input.

Cutting mAL_m8 alone:

| Map | Male P1 | Female | Meets primary? |
|---|---|---|---|
| A | 0,0,0 | WT | **No** |
| B | 11,13,12 | WT | Yes (and 3/3 male>0) |
| C | 1,0,1 | WT | Yes (2/3) |

So mAL_m8 is a robust **brake**, but not a robust **sufficient cut** under Map A.

---

## 2. Is AN09B017c’s gate identity only under c=−1?

**Yes — sign-dependent.**

| Map | sign | I_m | I_f | E_m | Role |
|---|---:|---:|---:|---:|---|
| A | −1 | 238 | 0 | 0 | **male-specific inhibitor** (rank #2) |
| B | +1 | 0 | 0 | **236** | **male-specific excitatory** (E only to male) |
| C | 0 | 0 | 0 | 0 | **silent** |

The “c is a last-hop inhibitory gate” claim exists **only on Map A**. On Map B
the same wiring is male-biased **excitation**. On Map C it does nothing.

---

## 3. Is there a small female-intact cut-set that does **not** depend on b/c sign?

**No (primary criterion on all three maps).**

Intersection of inhibitory types present in every map’s ranking includes
**mAL_m8, mAL_m2b, PVLP048, VES022** (always GABA/inhibitory). Testing that
ladder as a **sign-robust** cut-set:

| Robust set | Map A M | Map B M | Map C M | All-map primary? |
|---|---|---|---|---|
| {mAL_m8} | 0,0,0 | 11,13,12 | 1,0,1 | **No** (A fails) |
| {m8, m2b} | 0,0,0 | 15,19,15 | 2,3,1 | **No** |
| {m8, m2b, PVLP048} | 0,0,0 | 15,19,17 | 2,3,0 | **No** |
| {m8, m2b, PVLP048, VES022} | 0,0,0 | 18,16,17 | 1,6,2 | **No** |

**LABEL: none** — no set built only from types inhibitory in all maps
releases male on Map A while keeping female near WT.

### Per-map min-cuts (map-specific ranking)

| Map | Min cut (primary) | Types | Male | Female |
|---|---|---|---|---|
| A | top2 | **mAL_m8 + AN09B017c** | 0,1,1 | WT |
| A (secondary) | top3 | m8 + **c** + m2b | **5,5,1** | WT |
| B | top1 | mAL_m8 | 11,13,12 | WT |
| C | top1 | mAL_m8 | 1,0,1 | WT |
| C (secondary) | top4 | m8+PVLP048+m2b+VES022 | 1,6,2 | WT |

Map A **requires c**. That is a sign-dependent mechanism.

---

## 4. Architecture-level conclusions

| Statement | Kind |
|---|---|
| mAL_m8 delivers male-only I to P1 on **all** sign maps | **Sign-robust brake** |
| mAL_m2b, VES022, PVLP048 also male-biased I across maps | Sign-robust (weaker) |
| AN09B017c as inhibitory gate | **Sign-dependent** (only c=−1) |
| AN09B017c as excitatory male drive | Map B only |
| Female-intact male release without touching c (or flipping c) | **Fails on Map A** |
| Male expression is easy if c=+1 (WT already male-biased) | Sign sensitivity |
| Sex selectivity = few gates? | Partial: mAL_m8 is a real shared gate; **c’s role is not robust** |

### Sign-robust vs sign-dependent

```text
SIGN-ROBUST
  mAL_m8  ── male-specific delivered I on A,B,C
  mAL_m2b, VES022, PVLP048 ── male-biased I on A,B,C
  FLA001m/003m ── female-biased E on A,B,C

SIGN-DEPENDENT
  AN09B017c ── inhibitor only if −1; excitatory if +1; dead if 0
  Best Map-A cut-set (m8+c) ── requires c=−1
  Phenotype C without cut ── only on Map B
```

---

## 5. Core answers

1. **mAL_m8** — yes, male-specific inhibitory brake on all three maps.
2. **AN09B017c gate** — **only under c=−1**; not a sign-robust mechanism.
3. **Small female-intact cut-set independent of b/c sign** — **does not exist**
   under the primary criterion (Map A always needs c).
4. **Minimal stable combination** — none across maps. Closest robust brake
   is **mAL_m8 alone**, sufficient on B and C, **insufficient on A**.
5. **Sensitivity** — high: flipping b/c to +1 releases male P1 at WT
   (14,11,11); zeroing b/c nearly does (0,0,1). The package’s female-biased
   baseline and its “c as gate” story are **artifacts of the Glu→−1 map**.

---

## 6. Tracked types (delivered flux, male / female)

| Type | A | B | C |
|---|---|---|---|
| mAL_m8 I | 582 / 0 | 533 / 0 | 579 / 0 |
| mAL_m2b I | 246 / 55 | 304 / 55 | 209 / 55 |
| AN09B017c | I 238/0 | **E 236/0** | silent |
| VES022 I | 153 / 0 | 86 / 0 | 81 / 0 |
| PVLP048 I | 135 / 0 | 83 / 0 | 136 / 0 |
| FLA001m E | 1020/1330 | 1092/1330 | 1019/1330 |
| FLA003m E | 602/1009 | 568/1009 | 608/1009 |

---

## 7. Claim discipline

| Allowed | Not allowed |
|---|---|
| mAL_m8 is sign-robust male-specific I | mAL_m8 alone always releases male (false on A) |
| c gate is sign-dependent | c is an established inhibitory neuron in vivo |
| No all-map female-intact cut-set without c | “Male release is impossible” (Map B WT already male-biased) |
| Preference score / behavior | Mate preference claims |

## 8. Reproduce

```sh
uv run python run_sign_robust.py .experiment-data --out sign-robust-results.json
```
