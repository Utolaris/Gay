# Mechanism-search audit and next-experiment proposal

**Status:** analysis only. No new outcome-seeking runs. Pre-registered baselines
(`PROTOCOL.md`, `baseline-results.json`) and the post-hoc search
(`EXPLORATORY.md`) stand unchanged.

**Scope of prior claim (must not over-read):** Under the current MaleCNS graph +
bounded LIF + uniform candidate sensory encoding + eight-cell pC1_4a/b readout,
male > female was obtained only by neutralizing the female pathway. That is a
property of this model stack, not a claim about biological mechanism uniqueness.

---

## 1. Model incompleteness audit

| Layer | What the model has | What is missing / distorted |
|---|---|---|
| Connectome | Full classified v1.0 graph, anatomical synapse counts | No receptor identity, no synapse-level sign, no developmental state |
| Fast signs | ACh +1, GABA −1, Glu −1, other 0; **vAB3e/f forced +1**; **candidate_female forced +1** | Female cells are 13 GABA + 14 unclear at baseline; the override is what opens the female→P1 path |
| Inputs | Uniform Bernoulli Poisson, total 5550 Hz over all cells of a class | No labeled lines, no pheromone tuning, no per-cell identity, no adaptation, no concentration |
| Dynamics | Homogeneous LIF, one reversal, 300 ms | No neuromodulation (Fru/Pkd21 etc.), no persistent internal state, no arousal, no learning, no biomechanics |
| Readout | Eight literature-mapped pC1_4a/b | Not the genetic driver; direct vAB3→P1 wiring is only **8 synapses** in this graph (see `experiment/followup/MAPPING.md`) |
| Behavior | Spike-count preference score | No motor program, no choice, no persistence, no same-sex vs opposite-sex courtship distinction |

**Consequences for interpretation**

1. “Female-biased P1 response” is largely a **sign-map assumption**, not an
   emergent prediction from anatomy alone.
2. “Male disinhibition via mAL” is supported as a **conditional model result**
   (and literature-consistent with Kallman), but male absolute drive to P1 is
   anatomically weak.
3. Preference score ≠ mate preference. It cannot separate:
   - increased same-sex courtship,
   - bidirectional courtship,
   - true male-biased partner preference.

---

## 2. Signed pathway decomposition: candidate_male → P1

Data: prepared graph, baseline sign overrides (vAB3 + female excitatory).
Products are **anatomical count products**, not physiological efficacies.

### 2.1 No direct edge

- candidate_male → P1: **0 edges**
- candidate_female → P1: **0 edges**
- candidate_male → mAL: **0 edges**
- candidate_female → mAL: **0 edges**

Both cues must go through intermediates. mAL is not a direct sensory target.

### 2.2 Male 1-hop is local sensory, then the AN09B017 (vAB3-family) hub

Top male targets: LgLG6/7 (within-class), AN05B035, IN05B011a/b, then
AN09B017a–g. All male cells are ACh (+1).

### 2.3 Depth-2 male → intermediate → P1 is **dominated by inhibition**

| Intermediate | Type role | j_sign | path_sign | anatomical product |
|---|---|---|---|---|
| AN09B017c | Glu, not in vAB3e/f override | −1 | **−1** | 19329 |
| AN09B017b | Glu | −1 | **−1** | 3457 |
| AN09B017d | Glu | −1 | **−1** | 1940 |
| AN09B017e | vAB3e/f, forced +1 | +1 | **+1** | 1035 |
| AN03A008 | ACh | +1 | **+1** | 904 |
| AN09B017f | vAB3e/f, forced +1 | +1 | **+1** | 705 |
| AN09B017g | Glu | −1 | **−1** | 554 |

**Key structural fact:** the strongest male→P1 disynaptic anatomical routes are
**inhibitory glutamatergic AN09B017b/c/d**, not the canonical vAB3e/f excitatory
pair. That alone can keep male P1 sparse even when mAL output is blocked.

### 2.4 mAL is the largest but not the only P1 brake

Direct mAL → P1: **375 edges, 4251 synapses**, all sign −1.

Other strong inhibitory P1 sources (subset): SIP100m, SIP112m, SIP113m, VES022,
LH004m, PVLP048, GNG700m. Global mAL silence does not remove these.

Excitatory P1 sources include FLA001m, SIP105m, SIP025, pC1_1a (recurrent),
AVLP720m/721m/732m, VES206m, AN08B020 — shared central relays, not male-specific.

### 2.5 vAB3 story in this graph

- vAB3e/f → mAL: **5437 synapses** (very strong)
- vAB3e/f → P1 direct: **8 synapses** (almost none)
- Depth-2 vAB3→P1 is via mAL (inhibitory) and FLA001m / AVLP / pC1_3b (excitatory)

So “male cue → vAB3 → P1” is **not** a strong direct line in v1.0 anatomy. The
published functional pathway and this graph’s direct wiring disagree; any
mechanism claim must say which it is using.

### 2.6 Female override is load-bearing

- Female base signs: 14×0 (unclear), 13×−1 (GABA); **0 excitatory**
- After override: 27×+1
- Without override, female has no fast output → female P1 collapses (already
  seen in exploratory `no_fem_exc_*`)

Female 1-hop also converges on AN09B017g (strong) and shared AN05B035 / IN05B011
— so female and male share hubs; they are not parallel labeled lines in this
encoding.

---

## 3. Three phenotypes (do not conflate)

| Phenotype | Operational definition in this model | Current status |
|---|---|---|
| **A. Same-sex courtship increase** | Male-cue P1 (or downstream) rises under intervention while female response is irrelevant or also rises | **Supported in-model:** mAL output silence turns male P1 from 0 → 1–5 spikes (all seeds). Aligns with Kallman-style “mAL brake” story. |
| **B. Bidirectional courtship** | Both male and female cues drive readout above baseline under the same intervention | **Supported in-model:** mAL silence gives male 1–5 and female 12–22 simultaneously. |
| **C. True male-biased preference** | Male readout **exceeds** female under matched encoding, without destroying the female pathway | **Not supported** by pre-registered grid. Post-hoc “wins” required killing/neutralizing female output — that is pathway deletion, not preference. |

A behavioral layer (mounting, chasing, persistence) is absent; none of A–C is
a whole-animal courtship claim.

---

## 4. Routes that keep the female pathway intact

Ranked by **expected information gain** (what we learn), not by chance of
forcing a preferred score.

### E1 — Lesion map of P1 brakes under male vs female drive *(highest priority)*

- For each inhibitory class with substantial P1 weight (mAL subtypes, SIP100m,
  SIP112m, VES022, LH004m, …), apply `output_silence` **one class at a time**.
- Male and female cues; matched seeds; female pathway fully intact.
- Read out: ΔP1_male, ΔP1_female, global activity.
- **Question answered:** which non-mAL brakes actually gate male vs shared
  drive? Is mAL special, or one of several interchangeable brakes?

### E2 — Male-specific excitatory relay gain/activation

- Targets from decomposition with path_sign +1: **AN09B017e/f, AN03A008,
  AN09B002**.
- Interventions: `output_gain` ≤1 is only dampening; use `activate_tonic` /
  `activate_poisson` on these classes (not P1), or allow a documented
  `max_gain>1` sensitivity **labeled as amplification**.
- Keep female input and signs unchanged.
- **Question:** can a weak male excitatory branch be amplified enough to
  approach female drive without touching female?

### E3 — Selective relief of AN09B017b/c/d inhibition

- These are the strongest male→P1 inhibitory intermediates.
- Silence or reduce gain on **only those types** (not all Glu, not female).
- **Question:** is male sparseness dominated by this local glutamate brake
  rather than by mAL?

### E4 — Central E/I routing reweight on shared hubs

- Shared targets: AN05B035, IN05B011a/b. Both sexes hit them.
- Test whether shifting their outgoing E/I balance (bounded) changes male and
  female **differentially**.
- Risk: shared hub changes will move both; informative even if null.

### E5 — Recurrent / arousal state (persistent gain)

- Add a slow global or class-gated conductance (state variable), not a new
  sensory rate.
- **Question:** does a nonlinear state amplify sparse male drive more than the
  already-strong female drive? Likely amplifies both — measure ratio, not just
  absolute male spikes.

### E6 — Readout robustness

- Repeat key contrasts on **broad 49-cell pC1** and on literature P1-related
  subsets; do not stop at pC1_4a/b.
- **Question:** is “female bias” a readout-definition artifact?

### E7 — Input encoding audit (before more interventions)

- Sparse / cell-specific patterns inside LgLG6/7 vs LgLG5/8; matched total
  event count.
- **Question:** does uniform Poisson erase a male-labeled line that anatomy
  could support?

### E8 — Structural rewiring sensitivity *(last, clearly hypothetical)*

- Increase AN09B017e/f→P1 anatomical weights from 8 synapses toward
  literature-implied strength; keep everything else.
- **Question:** how much male drive is missing purely from this mapping
  discrepancy? Not a biological result — a “what if the paper’s wiring were
  in the graph” probe.

---

## 5. Proposed next round (concrete, pre-registerable)

Do **not** start by tuning preference. Run E1 then E2/E3 as a fixed grid.

### Protocol sketch (to freeze before running)

- Graph, LIF constants, sign map **unchanged** from baseline (including female
  excitatory override).
- Seeds 11, 12, 13; reversal −70 mV; matched male/female schedules.
- **Block 1 (E1):** WT vs silence of each of ~6–8 named inhibitory P1 sources;
  report full Δ tables.
- **Block 2 (E2/E3):** top 2 male-specific excitatory relays and top 2 male
  inhibitory intermediates from E1/decomposition; fixed intervention levels
  chosen **before** seeing P1 outcomes (e.g. tonic 6 mV, gain 0.25 / silence).
- Primary metrics: P1 male/female spikes, preference score, global spikes.
- **Explicit non-goals:** no direct P1 drive; no seed picking; no global
  threshold change; no claiming phenotype C unless male>female **with female
  pathway intact** on all seeds.
- Phenotype labels applied post hoc to whatever pattern appears (A / B / C /
  null).

---

## 6. Claim discipline

| Statement | Allowed? |
|---|---|
| In this model, mAL output silence unmasks a small male-cue P1 response | Yes (already measured) |
| In this model, male>female required neutralizing female fast output | Yes (exploratory) |
| mAL is the unique biological brake on male courtship | No — graph shows other P1 inhibitors; literature is broader |
| vAB3 is the main male excitatory relay to P1 in v1.0 anatomy | No — direct vAB3→P1 is 8 synapses; strongest male depth-2 paths are inhibitory AN09B017 |
| Any preference score here is mate preference | No |

---

## 7. Artifacts

- Pathway dump: run `uv run python analyze_pathways.py` (static, read-only).
- Baseline: `baseline-results.json` / `RESULTS.md`
- Post-hoc flip conditions: `EXPLORATORY.md` / `exploratory-male-bias.json`
