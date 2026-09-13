# AN09B017b/c sign + LgLG hop-3 routing (S-block)

**Scope:** MaleCNS v1.0 annotations, signed-path decomposition, and a frozen
female-intact intervention/sign-sensitivity grid. Not a parameter scan and not
a preference-tuning campaign.

Scripts: `analyze_bc_sign_hop3.py`, `run_bc_hop3.py`.
Protocol: [BC_HOP3_PROTOCOL.md](BC_HOP3_PROTOCOL.md).
Outputs: `bc-hop3-evidence.json`, `bc-hop3-results.json` (120/120 runs, 0 errors).
Seeds 11–13. Tests: 11 passed.

---

## 1. Is AN09B017b/c an inhibitory brake?

### 1.1 Annotation / neurotransmitter evidence

| Type | n | dimorphism | synonyms | consensus_nt | predicted_nt (conf) | receptorType | model sign |
|---|---:|---|---|---|---|---|---:|
| AN09B017a | 2 | male-specific | — | glutamate | glutamate (0.81) | empty | −1 |
| **AN09B017b** | 2 | male-specific | — | **glutamate** | **acetylcholine (0.94–0.96)** | **empty** | **−1** |
| **AN09B017c** | 2 | male-specific | — | **glutamate** | **acetylcholine (0.94)** | **empty** | **−1** |
| AN09B017d | 2 | male-specific | — | glutamate | acetylcholine (0.93–0.95) | empty | −1 |
| AN09B017e | 2 | male-specific | Yu 2010: vAB3 | glutamate | acetylcholine (0.95) | empty | −1 → **+1 override** |
| AN09B017f | 2 | male-specific | Yu 2010: vAB3 | glutamate | acetylcholine (0.94–0.96) | empty | −1 → **+1 override** |
| AN09B017g | 2 | male-specific | Yu/von Phillipsborn: vAB3 | glutamate | acetylcholine (0.92–0.97) | empty | −1 (no override) |

All AN09B017 subtypes are `fru_high` ascending neurons. Direct b/c→P1
synapses are small (b: 7, c: 55); b mainly targets mAL / AN05B035 / SIP, not P1.

### 1.2 Claim classes

| Claim | Kind |
|---|---|
| b/c are male-specific fru ANs; consensus NT is glutamate | **Anatomical fact** |
| LgLG6→b/c synapses ≫ LgLG5 (prior hop-2) | **Anatomical fact** |
| receptorType is absent for every AN09B017 cell | **Anatomical fact (missing data)** |
| Fast sign glutamate → −1 | **Model assumption** |
| e/f forced +1 because literature vAB3 | **Model assumption** (literature-anchored) |
| g has vAB3 synonyms but stays −1 | **Model inconsistency** |
| predicted_nt = ACh for b/c | **Competing prediction, unused by model** |
| b/c inhibit P1-related targets in vivo | **Biological hypothesis — not established** |
| postsynaptic GluCl vs AMPA-like receptor at b/c→P1 sources | **Unknown in this dataset** |

**Verdict:** calling b/c an “inhibitory brake” is a **model assumption layered
on consensus glutamate**, not a receptor-validated anatomical fact. The
MaleCNS NT-prediction column actively disagrees (ACh). No postsynaptic
receptor annotation exists here to break the tie.

### 1.3 Anatomical sign sensitivity (static)

Hop-2 weighted P1-source access (net E−I), by b/c sign:

| Variant | LgLG6 net | LgLG5 net |
|---|---:|---:|
| b/c = −1 (baseline) | **−21441** | −1583 |
| b/c = +1 | **+20075** | −1393 |
| b/c = 0 | −683 | −1488 |

Hop-3 signed product mass to P1:

| Variant | LgLG6 net | LgLG5 net | ranking |
|---|---:|---:|---|
| b/c = −1 | +3.88e7 | **+7.98e7** | female ~2× |
| b/c = +1 | **+1.26e8** | +7.06e6 | **male ~1.8×** |
| b/c = 0 | +8.15e7 | +7.52e7 | near-tie |

**Flipping b/c to +1 reverses the entire hop-2/3 drive ranking.** The
“female-advantage pathway mass” is conditional on Glu→−1 for b/c.

### 1.4 LIF sign sensitivity (S1, female intact, 8×2×3=48 runs)

| Condition | Male P1 | Female P1 | pref | C? |
|---|---|---|---:|---|
| base_WT | 0,0,0 | 11,8,8 | −1.00 | no |
| base_mAL | 4,5,1 | 18,18,20 | −0.70 | no |
| **bc_exc_WT** | **14,11,11** | **11,8,8** | **+0.15** | **YES** |
| bc_exc_mAL | 30,28,20 | 18,18,20 | +0.16 | no (seed13 tie 20=20) |
| bc_zero_WT | 0,0,1 | 11,8,8 | −0.93 | no |
| bc_zero_mAL | 10,11,3 | 18,18,20 | −0.42 | no |
| **bcg_exc_WT** | 14,11,11 | **0,0,0** | +1.00 | **YES** |
| **bcg_exc_mAL** | 30,28,20 | 11,16,15 | +0.29 | **YES** |

**Key result:** under baseline signs the model is female-biased. Treating
b/c as **excitatory** (the predicted_nt direction) is enough to produce
phenotype C with the female pathway fully intact and no other change
(`bc_exc_WT`: male 14,11,11 > female 11,8,8 on all seeds).

Flipping **g** as well (`bcg_exc`) collapses female to 0 at WT: g is heavily
loaded by LgLG5 and projects mostly onto mAL, so excitatory g would amplify
the shared brake rather than the female→P1 drive. g’s literature vAB3 synonym
does **not** make “force all AN09B017 +1” a safe reading.

---

## 2. Why LgLG5 hop-3 positive mass ≈ 2× LgLG6

Depth-3 paths: source → t1 → t2 → P1. Baseline signs.

### 2.1 Mass budget

| | LgLG6 | LgLG5 |
|---|---:|---:|
| n paths | 198548 | 132274 |
| positive mass | 3.50e8 | 1.93e8 |
| negative mass | −3.11e8 | −1.13e8 |
| **net** | **+3.88e7** | **+7.98e7** |

Female does not win by having more positive raw mass — male has more of
both. Female wins because male’s **extra negative mass** is larger than
male’s extra positive mass.

### 2.2 Ranked roles (by hop-3 mass as t1)

**Female-side amplifier / male-specific sink**

| t1 type | LgLG6 mass | LgLG5 mass | Δ(5−6) | role |
|---|---:|---:|---:|---|
| **AN09B017c** | **−4.68e7** | −2.1e5 | +4.66e7 | **male-specific inhibitory sink** |
| **AN09B017g** | +5.58e6 | **+3.80e7** | +3.24e7 | **female-biased positive carrier** |
| LgLG6 (recurrence) | −5.11e6 | −9.5e3 | +5.10e6 | male within-class negative loop |
| AN03A008 | −1.81e6 | −4.6e4 | +1.76e6 | male-loaded, net negative at depth 3 |

**Male-heavier mass (not all “sinks”)**

| t1 type | LgLG6 | LgLG5 | note |
|---|---:|---:|---|
| AN09B017d | **+5.63e7** | +3.74e7 | shared, **positive** product — disinhibitory motif under −1 sign |
| AN09B017a | +2.82e7 | +1.81e7 | shared positive |
| AN05B035 | +1.80e7 | +1.18e7 | shared hub, positive |
| AN09B017f | −1.69e7 | −2.21e7 | forced +1; female more negative |
| IN05B011a | +2.44e6 | +1.2e5 | male-specific load |

### 2.3 Minimal central routing switch

The first irreversible divergence is not a shared hub. It is the
**LgLG6→AN09B017c labeled-line edge** (≈170× per-cell vs LgLG5), which
injects the dominant negative hop-3 mass in male. LgLG5 instead loads
**AN09B017g**, whose depth-3 product is strongly positive.

There is no single shared relay whose E/I reweight separates the sexes
without also moving both (E4/S2 hub tests). The switch is **edge-selective
routing into c vs g**, not a global gain on AN05B035 / IN05B011.

---

## 3. Female-intact local interventions (S2, 72 runs)

References: WT male 0,0,0 / female 11,8,8; mAL male 4,5,1 / female 18,18,20.

| Condition | Male P1 | Female P1 | pref | female == mAL-ref? | Δmale vs mAL |
|---|---|---|---:|---|---|
| mAL | 4,5,1 | 18,18,20 | −0.70 | — | — |
| **c_mAL** | **10,6,3** | 18,18,20 | −0.51 | **yes** | +6,+1,+2 |
| **bc_mAL** | **10,11,3** | 18,18,20 | −0.42 | **yes** | +6,+6,+2 |
| bcd_mAL | 9,9,5 | 18,18,20 | −0.42 | yes | +5,+4,+4 |
| a008_mAL | 1,2,1 | 18,18,20 | −0.87 | yes | −3,−3,0 |
| a008_bc_mAL | 8,6,4 | 18,18,20 | −0.52 | yes | +4,+1,+3 |
| a008_c_mAL | 7,10,4 | 18,18,20 | −0.46 | yes | +3,+5,+3 |
| in11a_mAL | 1,1,0 | 18,19,20 | −0.93 | no (F↑) | −3,−4,−1 |
| in11a_bc_mAL | 10,9,5 | 18,19,20 | −0.41 | no (F↑) | +6,+4,+4 |
| d_g25_mAL | 0,3,0 | 16,18,20 | −0.91 | no (F↓) | **−4,−2,−1** |
| a008_in11a_mAL | 0,0,0 | 18,19,20 | −1.00 | no (F↑) | −4,−5,−1 |

### Findings

1. **Minimal local mechanism (frozen signs):** `mAL_silence + AN09B017c
   silence` (optionally + b). Male roughly doubles-to-triples vs mAL alone;
   **female is bit-identical** to the mAL-alone female response. This matches
   the hop-2 anatomy: c is the male-specific sink.
2. **Including d is not better** (bcd_mAL ≤ bc_mAL). d carries large
   *positive* hop-3 mass under the current sign map; partial reduction
   (`d_g25`) **hurts** male. d is not a brake to peel off.
3. **AN03A008 tonic 6 mV does not add** to bc/c relief (a008_bc_mAL ≤
   bc_mAL). Male excitatory relay amplification is not the bottleneck once
   mAL and c are removed — or 6 mV is too weak to matter at this readout.
4. **IN05B011a is not a male gate** in this encoding: silence alone is null
   or harmful; combined with bc it is no better than bc alone.
5. **No local frozen-sign condition reaches phenotype C.** Best male intact
   drive remains 10,11,3 vs female 18,18,20 (pref ≈ −0.42).

---

## 4. Does a real central preference-reversal circuit exist?

Answer, separated by evidence class:

| Statement | Status |
|---|---|
| Connectome routes LgLG6 into AN09B017c ≫ LgLG5 | **Yes (anatomy)** |
| That routing can gate male P1 when mAL is off | **Yes (model, S2)** |
| That routing **reverses** male vs female under baseline signs | **No** |
| Reversal appears if b/c are excitatory (predicted ACh) | **Yes (S1, model assumption)** |
| b/c are excitatory in vivo | **Not established** |
| A minimal local intervention under frozen signs flips preference | **No** |
| Preference score here is behavioral mate choice | **No** |

**Synthesis.** A central **routing asymmetry** exists and is real anatomy:
LgLG6→c vs LgLG5→g. Under the package’s frozen Glu→−1 map it is a male
brake, not a male-biased drive line — releasing it (with mAL) raises male
without touching female, but female remains ~2× higher. A true preference
reversal in this model stack appears only when the **sign of b/c is changed**
to the predicted-ACh direction, which is a hypothesis about receptors, not a
demonstrated circuit fact. Without receptor-level evidence, phenotype C from
`bc_exc` is a **sensitivity result**, not a mechanism claim.

Structural what-if E8 (vAB3e/f→P1 ×100) remains the other male-approach
path, and is explicitly a mapping-discrepancy probe.

---

## 5. What was not done (and why)

- No expansion of tonic/gain ladders after seeing null a008.
- No silencing of g or of the female sensory class.
- No global threshold / weight retuning.
- No seed selection; all 120 runs reported.
- No postsynaptic receptor data inventing — field is empty in MaleCNS v1.0
  tables used here.

---

## 6. Reproduce

```sh
uv run python analyze_bc_sign_hop3.py
uv run python run_bc_hop3.py .experiment-data --workers 6
uv run pytest -q
```

## 7. Next experimental need (not a parameter scan)

To close the sign question without more model grids: **postsynaptic receptor
identity at AN09B017b/c → P1-source synapses** (GluCl vs AMPA/kainate-like),
or functional paired recordings. Until then, both “b/c brake” and
“b/c ACh relay” must be reported as competing model maps.
