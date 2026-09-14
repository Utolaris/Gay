# P1_CURRENT autopsy results

Protocol: [P1_CURRENT_PROTOCOL.md](P1_CURRENT_PROTOCOL.md).
Raw: `p1-current-results.json`. Seeds 11–13. Female path intact.
Attribution uses **delivered weighted flux onto P1**, not source spike counts.

WT P1: male **0,0,0** · female **11,8,8**.

---

## 1. What actually hits P1 (Phase A)

Per-P1-cell traces (Vm, ge, hi, E/I flux, first-spike latency) are in
`phase_a_full_seed11`. Male never reaches threshold; female spikes on all
8 cells over the window.

### Delivered inhibition onto P1 (mean weight-flux units)

| Type | I_male | I_female | ΔI | male-specific score |
|---|---:|---:|---:|---:|
| **mAL_m8** | **582** | **0** | +582 | **582** |
| **AN09B017c** | **238** | **0** | +238 | **238** |
| **mAL_m2b** | **246** | 55 | +191 | 157 |
| VES022 | 153 | 0 | +153 | 153 |
| PVLP048 | 135 | 0 | +135 | 136 |
| SIP100m | 94 | 16 | +78 | 65 |
| AVLP712m | 57 | 0 | +57 | 57 |
| LH004m | 61 | 12 | +49 | 41 |
| SIP113m | 230 | 342 | −112 | −45 |
| SMP165 | 133 | 378 | −245 | −64 |
| mAL_m9 | 558 | 741 | −183 | −79 |

**mAL_m8 and AN09B017c deliver large I to P1 only under male cue** (I_f = 0).
They are last-hop **male-specific** inhibitory inputs in the delivered-current
sense (not merely “they spike more”).

Female-heavier inhibitors (mAL_m9, SMP165, SIP113m) are shared/female-side
brakes — cutting them would not selectively help male.

### Delivered excitation onto P1

| Type | E_male | E_female | ΔE |
|---|---:|---:|---:|
| FLA001m | 1020 | 1330 | **−310** |
| FLA003m | 602 | 1009 | **−406** |
| pC1_3b | 217 | 329 | −112 |
| LH006m | 169 | 268 | −99 |
| pC1_3c | 90 | 21 | +69 |

Female advantage is **both** less male-specific I **and** stronger shared E
(FLA001m/FLA003m). Sex selectivity is not only an inhibitory story.

---

## 2. Minimal causal cut-set (Phase B)

Ablation = `gain=0` on listed types (sources may still spike; delivered
current is blocked). Ranking frozen from Phase A scores.

| Condition | Types | Male P1 | Female P1 | male>0 seeds | female ≈ WT |
|---|---|---|---|---:|---|
| top1 | mAL_m8 | 0,0,0 | 11,8,8 | 0 | 3/3 |
| **top2** | **mAL_m8 + AN09B017c** | **0,1,1** | **11,8,8** | **2/3** | **3/3** |
| top3 | + mAL_m2b | 5,5,1 | 11,8,8 | 3/3 | 3/3 |
| top2 + full mAL | m8+c + mAL class | 10,6,3 | 18,18,20 | 3/3 | 3/3 (but F↑) |

**Pre-registered min-cut: `top2`** (male P1>0 on ≥2/3 seeds, female within
WT−2 on ≥2/3 seeds).

`top3` is the robust all-seeds male release with **bit-identical female**.
Adding the rest of mAL raises female — not male-specific.

---

## 3. Answers

### Q1 — Which inhibitory currents hold male P1 below threshold?

Delivered-current ranking: **mAL_m8** (582 units, male-only) and
**AN09B017c** (238, male-only), then **mAL_m2b**. These are the last-hop
I that male cue actually injects into P1 and female does not.

### Q2 — Male-specific vs shared brakes?

| Class | Examples |
|---|---|
| **Male-specific I** | mAL_m8, AN09B017c, VES022, PVLP048 (I_f=0) |
| Shared / female-heavier I | mAL_m9, SMP165, SIP113m |
| Female-biased E | FLA001m, FLA003m |

### Q3 — Small inhibitory cut-set that releases male, spares female?

**Yes.** Two types (`mAL_m8` + `AN09B017c`) give male P1 0→1 on 2/3 seeds
with female **exactly WT**. Three types give male 5,5,1 with female still
WT. This is the smallest demonstrated female-intact male release under
unchanged signs.

### Q4 — Distributed code vs gate neurons?

**Both, at different scales.**

- A **2–3 type cut-set** already unmasks male spikes (gate-like).
- Full latent drive still requires removing most P1 I (prior autopsy: 267
  spikes only when ~all P1 inhibition is cut) — remaining inhibition is
  distributed.
- Female E advantage (FLA*) is parallel and untouched by the male-specific
  cut, which is why female stays high.

So sex selectivity here is **not** a single gate neuron, but it is also
**not** a featureless distributed code: a short male-specific inhibitory
list sits at the last hop.

---

## 4. Method note

`output_silence` / gain=0 does **not** stop source spiking. Phase B measures
P1 spikes after blocking delivery. Phase A ranks by **arriving weighted
flux** (`i_by_type`), which is zero if delivery is blocked even when the
source still fires.

---

## 5. Claim discipline

| Statement | Status |
|---|---|
| mAL_m8 delivers I to P1 only under male cue (this encoding) | Model, delivered-current |
| AN09B017c is a last-hop male-specific I source to P1 | Model, delivered-current |
| top2/top3 release male P1 with female intact | Holds (seeds reported) |
| These types are the in vivo male gates | Not established |
| Preference / courtship behavior | Not claimed |

## 6. Reproduce

```sh
uv run python run_p1_current.py .experiment-data --out p1-current-results.json
```
