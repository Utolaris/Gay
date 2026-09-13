# AN09B017b/c output sign adjudication

**Question:** is the real synaptic output of AN09B017b/c excitatory or inhibitory?

**Short answer:** **cannot be settled to a single physiological fact from MaleCNS v1.0 tables alone.** Ranking the competing maps, the **literature-consistent circuit role of the AN09B017/vAB3 class is excitatory**; the package’s Glu→−1 “brake” for b/c is a **weaker default** than treating them like e/f. Subtype **a** is the only AN09B017 cell whose EM NT predictor also says glutamate.

Local sources: `nt.feather`, `annotations.feather`. Prior sensitivity: `bc-hop3-evidence.json`, S1 in `bc-hop3-results.json`.

---

## 1. Transmitter identity (what is released)

| Layer | AN09B017b/c | Note |
|---|---|---|
| `consensus_nt` | **glutamate** | |
| `ground_truth` | **glutamate** | identical to consensus for all 85,484 cells that have it (0 disagreements) |
| `predicted_nt` | **acetylcholine**, conf 0.94–0.96 | body-level EM predictor |
| `celltype_predicted_nt` | **acetylcholine**, conf ~0.95 | type-level |
| Literature synonym on b/c | **none** | e/f: “Yu 2010: vAB3”; g: also vAB3 |
| `mancType` | **AN09B017** | one MANC type for a–g |
| `hemibrainType` | empty | no hemibrain cross-map |
| `receptorType` | **empty** | only ppk23/25/IR52b exist package-wide — pheromone channels, not synaptic Glu/ACh receptors |

`ground_truth` is **not independent of consensus**: when present, consensus copies it. The real conflict is:

> consensus/ground_truth = Glu  vs  EM predictor = ACh

Across the whole MaleCNS NT table only **77 cells** are `gt=Glu` + `pred=ACh`. AN09B017 **b–g** (12 cells) are a large share of that rare disagreement set. Subtype **a** (2 cells) is the only AN09B017 member where predictor also says Glu.

---

## 2. Sign (excitatory vs inhibitory) — what is missing

Even if transmitter is Glu, **sign is not determined**. Stürner et al. (Nature 2025) state the gap explicitly:

> “Neurotransmitter predictions … lack neuropeptide predictions and receptor expression data, an important gap given that neurotransmitters such as glutamate can be excitatory or inhibitory.”

There is **no postsynaptic receptor annotation** for b/c→P1-source synapses in this dataset. GluCl-type ⇒ inhibitory; AMPA/kainate-like ⇒ excitatory. Unknown here.

If the ACh prediction is correct, sign is unambiguously +1 under the package’s ACh→+1 map — no Glu dual-sign problem.

---

## 3. Class identity: AN09B017 = vAB3

| Source | Statement |
|---|---|
| Stürner et al., Nature 2025 | “AN09B017 seems to match *fru*-expressing **vAB3** neurons, which transmit contact pheromone signals from the front legs to **P1**” |
| Rubin et al., Curr Biol 2026 | “male-specific ascending neurons, **AN09B017 (previously referred to as vAB3)** … responsive to 7,11-HD” |
| Ryba et al., Curr Biol 2026 | circuit diagram: “**excitatory vAB3** neurons and inhibitory mAL neurons to P1” |
| MaleCNS synonyms | only **e/f** carry “Yu 2010: vAB3”; **g** also vAB3; **b/c/d/a** have empty synonyms |

So literature treats **AN09B017 as a class** ≈ vAB3 and **vAB3 as excitatory** onto P1. The package forced only the four e/f body IDs to +1 because only those cells had Yu-2010 string synonyms — that override is **under-applied relative to class-level identity**.

---

## 4. Competing maps, ranked

| Rank | Map | Predicted sign of b/c | Support | Problem |
|---:|---|---:|---|---|
| 1 | **Class-level vAB3 / circuit literature** | **+1 excitatory** | AN09B017 = vAB3; “excitatory vAB3→P1”; same MANC type as e/f which model already forces +1 | no cell-resolved physiology on b vs e |
| 2 | **EM NT predictor (ACh)** | **+1** | conf ~0.95; whole b–g type; rare disagreement with Glu gt | model not trained only on these cells; no receptor check |
| 3 | **consensus Glu + GluCl default** | **−1 inhibitory** | consensus/gt Glu; package Glu→−1; E1/E3 causal gate under this map | contradicts ACh prediction; Glu dual-sign; no receptors; only 77-cell conflict class |
| 4 | Subtype-heterogeneous | a: −1; b–g: +1 | a alone predicts Glu; b–g predict ACh | untested; morphology a–g may still share NT |

---

## 5. What the model got wrong (or at least over-committed)

In `PROTOCOL.md` / `prepare_data.py`:

- Fast signs: ACh +1, GABA −1, **Glu −1**.
- Overrides: vAB3e/f bodies forced +1; female sensory forced +1.
- **b/c/d/g stay −1** despite being the same MANC type `AN09B017` and despite g carrying vAB3 synonyms.

Consequences already measured (S1):

- With b/c = −1 (package default): male hop-2 access −21441; LIF female-biased.
- With b/c = +1: hop-2 access +20075; **`bc_exc_WT` → phenotype C** (male 14,11,11 > female 11,8,8, female intact).

So the “male-specific inhibitory sink” story is **conditional on a sign choice that literature and the EM predictor both push against**.

---

## 6. Adjudication

| Statement | Status |
|---|---|
| b/c release glutamate | **Dataset consensus** (not receptor-proven) |
| b/c release acetylcholine | **EM predictor, high confidence** — competing claim |
| b/c are inhibitory onto P1-related targets | **Not established**; model assumption |
| b/c are excitatory like other vAB3 | **Most consistent with class identity + circuit literature** |
| Postsynaptic receptor at b/c→P1 sources | **Unknown in this dataset** |
| Package’s “bc brake” is a biological fact | **No** |
| Package’s “bc brake” is a valid model sensitivity | **Yes**, and it should be labeled as such |

### Most likely answer (with residual uncertainty)

**Treat AN09B017b/c as excitatory (+1) for circuit interpretation**, same as e/f, unless receptor-level data show GluCl at their major postsynaptic targets.

Confidence: **moderate** for excitatory class role; **low** for “definitely ACh vs Glu-with-AMPA”; **high** that the current −1 brake is **not** a stronger fact than +1.

Subtype caveat: **AN09B017a** remains the best Glu candidate within the class (predictor agrees with consensus). If any AN09B017 subtype is a true inhibitory Glu neuron, it is more likely **a**, not b/c.

---

## 7. What would settle it (not another parameter scan)

1. **Postsynaptic receptor identity** at AN09B017b/c→(mAL, SIP105m, AVLP, P1 sources): GluClα vs GluR vs nAChR (transcriptomics / immunostaining / connectome-constrained receptor map if released).
2. **Functional**: optogenetic activation of b/c (split-GAL4 if available) + P1 or mAL imaging; sign of the PSP.
3. **Vesicle / NT immuno** specifically on AN09B017b/c terminals, not class-level consensus.
4. Until then, every result that depends on b/c sign must be reported under **both maps** (already done in S1).

---

## 8. Claim discipline for this package

| Allowed | Not allowed |
|---|---|
| “Under Glu→−1, b/c act as a male hop-2 brake in this LIF model” | “AN09B017b/c inhibit P1 in vivo” |
| “Flipping b/c to +1 produces phenotype C with female intact” | “Male-biased courtship is explained by b/c ACh release” |
| “AN09B017 class matches vAB3; literature calls vAB3 excitatory” | “Only e/f are vAB3; b/c are a different cell type” |
| “Receptor data are missing; Glu dual-sign is unresolved for consensus-Glu cells” | Choosing −1 because it fits a female-biased baseline |

---

## 9. Reproduce (local tables)

```sh
uv run python - <<'PY'
# see bc-hop3-evidence.json evidence_ledger and claim_classes
import json
print(json.load(open('bc-hop3-evidence.json'))['claim_classes'].keys())
PY
```

Literature anchors (not in-repo):

- Stürner et al. 2025, Nature 643:158–172 (AN09B017 ≈ vAB3; Glu dual-sign gap).
- Rubin et al. 2026, Curr Biol (AN09B017 previously vAB3; 7,11-HD).
- Ryba et al. 2026, Curr Biol (excitatory vAB3 → P1; inhibitory mAL).
- Clowney et al. 2015, Neuron (vAB3 pheromone→P1 circuit; cited by Stürner).
