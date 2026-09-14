# Male sensory × reward/modulatory convergence onto courtship

Static MaleCNS v1.0 anatomy. Script: `analyze_convergence.py`.
Output: `convergence-report.json`.

**Short answer:** there is **no robust hop-1 site** where male pheromone GRN
drive and dopamine reward (PAM/PPL) meet on a courtship-projecting neuron.
Octopamine makes only token contacts on the first male central hub
(AN09B017). Real OA/PPL-rich courtship relays sit **downstream**, reached by
male cue only at **hop 2**, mostly via AN09B017 / AN03A008.

---

## 1. Populations

| Set | n |
|---|---:|
| Male LgLG6/7 | 37 |
| PAM (DA) | 316 |
| PPL (DA) | 24 |
| OA | 37 |
| Courtship core (pC1* ∪ mAL ∪ P1) | 231 |

## 2. Hop-1: male targets ∩ modulatory targets

| | |
|---|---:|
| Male hop-1 targets | 253 |
| Modulatory hop-1 targets | 74,943 |
| **Overlap** | **50** |

All 50 overlapping cells can reach P1 within ≤2 hops. They are almost all
**AN09B017\*** and **AN05B\*** ascending neurons — the same male sensory hub
already mapped in S-block.

### Modulatory weight onto that hub (per type, Σ synapses)

| Type | PAM | OA | PPL |
|---|---:|---:|---:|
| AN09B017a | 0 | 0 | 0 |
| AN09B017b | 0 | 4 | 0 |
| AN09B017c | 0 | 8 | 0 |
| AN09B017d | 0 | 13 | 0 |
| AN09B017e/f/g | 0 | 2–7 | 0 |
| AN05B035 | 0 | 2 | 0 |
| AN03A008 | 0 | 8 | 0 |
| IN05B011a | 0 | 0 | 0 |

**PAM never innervates the male sensory first-order hub.** OA does, but at
token scale (single-digit to low tens of synapses).

### Courtship core with hop-1 male **and** modulatory input

**n = 0.** No pC1 / mAL / P1 cell receives direct male GRN input together
with PAM/PPL/OA.

---

## 3. Hop-2: male → AN09B017/AN03A008 → mod-recipient courtship relays

Among 1,097 modulatory-recipient cells that are courtship-core or →P1,
**948** also have male hop-2 access. Top by male×mod×P1:

| Type | male hop-2 | PAM | OA | →P1 | Dominant male via |
|---|---:|---:|---:|---:|---|
| **VES022** | 5–10×10³ | 0 | **155–217** | 45–81 | **AN03A008** |
| SIP025 | 1.3×10⁵ | 0 | 10 | 216–281 | AN09B017c/b |
| SIP105m | 3.4×10⁵ | 0 | 1 | 285 | AN09B017b, AN05B035 |
| SIP106m | 4–5×10⁴ | 1–5 | 6–22 | 90–92 | AN09B017e/f/c |
| mAL_m1/m2b/m5* | 0.6–2×10⁵ | 0 | 2–7 | 66–176 | AN09B017b/a, AN05B035 |
| AVLP721m | 4–6×10⁴ | 0 | 11–14 | 133–153 | AN09B017c/e |
| pC1_11a | 3.5×10⁴ | 0 | 17 | 72 | AN09B017b/c/e |

Pure-modulatory courtship relays (heavy DA/OA, **zero** male hop-1):
SIP106m, SMP163, OA-VPM3, OA-VUMa6, AOTU100m, pC1_1a, FLA001m — these are
reward/arousal→courtship lines that do **not** sample male GRNs directly.

---

## 4. Where convergence actually is

```text
Male GRN (LgLG6/7)
    │ hop-1
    ▼
AN09B017b/c/d/e/f, AN05B035, AN03A008     ← PAM=0, OA token only
    │ hop-2
    ▼
VES022 / SIP025 / SIP105m / SIP106m / mAL / AVLP / pC1_*
    ▲
    │ hop-1
OA (strong on VES022, SIP106m, SMP163)    ← PAM still thin here
PPL / PAM                                  ← mostly other branches
```

Ranked **biologically plausible coincidence detectors** (not P1-direct):

1. **VES022** — strongest OA among male-reachable courtship relays; male via
   AN03A008 (excitatory). Best candidate for OA×pheromone integration.
2. **SIP106m** — moderate OA + male via vAB3e/f/c + direct →P1.
3. **SIP025 / SIP105m / AVLP721m** — male-dominant, light OA; sensory relays
   with weak modulatory bias, not reward hubs.
4. **mAL subtypes** — male hop-2 is huge (so male can *drive the brake*), OA
   is token; not a reward-learning site.
5. **AN09B017 itself** — first male central hub; **not** a DA target. This is
   why MB/PAM pairing never touched the pheromone stream in REWARD block.

---

## 5. Implications for reward learning

| Prior assumption | Anatomical status |
|---|---|
| Male cue and PAM meet on KCs | **False** — LgLG→KC = 0 |
| Male cue and PAM meet on AN09B017 | **False** — PAM=0 onto AN09B017 |
| Male cue and OA meet on AN09B017 | **Token only** (2–13 syn) |
| OA can modulate male→courtship at VES022 / SIP106m | **Supported** (hop-2) |
| DA reward can gate the male pheromone line at hop-1 | **Not supported** |

So “reward association on the natural male→courtship line” would have to
happen at **VES022 / SIP106m-like** second-order relays under octopamine,
not at PAM–KC or PAM–AN09B017. That is a different experiment from the
MB block and still would not bypass mAL / b/c signs for P1 expression.

---

## 6. Claim discipline

| Statement | Kind |
|---|---|
| No hop-1 male∩PAM cell on courtship core | Anatomical |
| AN09B017 receives no PAM, only token OA | Anatomical |
| VES022 is OA-rich and male-reachable via AN03A008 | Anatomical |
| VES022 is the in vivo OA×pheromone integrator | **Hypothesis** |
| This explains REWARD-block nulls | Consistent, not proof of causation |

## 7. Reproduce

```sh
uv run python analyze_convergence.py
```
