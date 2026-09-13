# Why LgLG5 effective drive > LgLG6: signed pathway decomposition

**Scope:** anatomical findings on MaleCNS v1.0 + bounded LIF sensitivities.
Not a preference-score campaign. Female pathway intact throughout.

Scripts: `analyze_lglg_pathways.py`, `analyze_lglg_bridge.py`.
Outputs: `pathway-decomp-lg6-lg5.json`, `pathway-bridge-lg6-lg5.json`,
`lglg-bottleneck-confirm.json`.

---

## 1. Anatomical findings

### 1.1 Populations

| Type | n | Fast sign (baseline) |
|---|---:|---|
| LgLG6 (male GRN) | 16 | all ACh +1 |
| LgLG5 (female GRN) | 13 | base Glu/unclear; **forced +1 in model** |
| LgLG7 / LgLG8 | 21 / 14 | same pattern (7: +1; 8: forced +1) |

Neither class projects **directly** to the 8-cell P1 readout (0 edges).

### 1.2 Hop-by-hop mass (LgLG6 vs LgLG5)

| Hop | LgLG6 net E−I syn | LgLG5 net E−I syn | Note |
|---:|---:|---:|---|
| 1 | +10003 (E only) | +6045 (E only) | similar fan-out 45 vs 40 |
| 2 | **−65607** (E/I=0.70) | **−85312** (E/I=0.56) | both net-inhibitory |
| 3 | +3.4e6 | +2.5e6 | frontier truncated |
| 4 | ~0 | ~0 | near-balanced global mix |

Raw hop-2 E/I does **not** explain female advantage by itself (LgLG5 is
*more* net-inhibitory at hop 2 in absolute synapse count). The divergence is
in **which P1-relevant classes** receive that mass.

### 1.3 First major divergence: hop-2 access to P1 source classes

Weighted anatomical access = Σ syn(source→T) × syn(T→P1), by T’s fast sign.

| | LgLG6 | LgLG5 | ratio 6/5 |
|---|---:|---:|---:|
| hop2 → **P1-excitatory** sources | 1649 | 1075 | 1.53 |
| hop2 → **P1-inhibitory** sources | **23090** | **2658** | **8.69** |
| **Net (E−I) P1-source access** | **−21441** | **−1583** | **13.6× more negative for male** |

**This is the first hop where effective drive to the P1 microcircuit diverges
strongly: hop 2, via preferential male loading onto inhibitory P1 sources.**

### 1.4 Male-specific bottleneck (anatomy)

Per-cell synapses from GRN class into type T, and T→P1 synapses:

| Type | sign | 6/cell →T | 5/cell →T | ratio | T→P1 syn |
|---|---|---:|---:|---:|---:|
| **AN09B017c** | −1 | **39.4** | 0.2 | **171×** | 55 |
| **AN09B017b** | −1 | **60.7** | 0.2 | **263×** | 7 |
| IN05B011a | −1 | **70.9** | 4.5 | **16×** | 0 (indirect brake hub) |
| AN05B035 | −1 | 92.4 | 71.8 | 1.3× | 0 |
| AN09B017d | −1 | 33.3 | 28.0 | 1.2× | 7 |

LgLG6 preferentially innervates **Glu AN09B017b/c** (and strongly IN05B011a).
LgLG5 almost does not. These are the male-specific hop-2 inhibitory sinks.

Signed 2-hop product source→T→P1:

| Via | LgLG6 prod | LgLG5 prod | ratio |
|---|---:|---:|---:|
| AN09B017c | −17328 | −83 | **209×** |
| AN09B017b | −3430 | −12 | **286×** |
| AN09B017d | −1940 | −1264 | 1.5× |
| AN09B017f (E) | +705 | +1005 | 0.70× |
| AN03A008 (E) | +904 | +24 | 38× |

Male is not short on *all* excitatory hop-2 (AN03A008/f are fine); it is
flooded into b/c inhibition.

### 1.5 Female-side structure

- LgLG5 → AN09B017**g** (−1): 22.2 syn/cell vs LgLG6 2.7 (8×); but g→P1 is only 9 syn.
- LgLG5 → AN09B017**f** (+1, vAB3): 53.9/cell vs 33.7 — modest excitatory edge.
- Neither class reaches mAL or FLA001m/SIP105m directly (0 syn).
- LgLG7 vs LgLG8: LgLG8 loads AN09B017g even harder (bridge 188M vs 14M);
  **LgLG8-only drive is silent in the model** (E7) — g-heavy ≠ P1 drive here.
  Female sensory efficacy in this encoding is carried by **LgLG5**, not LgLG8.

### 1.6 Shared hubs (both sexes)

AN05B035, IN05B011a/b, AN09B017a/d/e/f/g, AN05B023*, IN04B079.
Hop-1 target overlap LgLG6∩LgLG5 = 64 cells. No male-exclusive labeled line
into the courtship core at hop 1.

### 1.7 Recurrence

Both classes return to themselves within 3 hops (source_in_3hop_out = 1.0).
No meaningful recurrence asymmetry at this resolution.

### 1.8 Signed path products to P1 (anatomical, depth 1–4)

| Depth | LgLG6 signed product | LgLG5 | hits 6/5 |
|---:|---:|---:|---|
| 2 | −21441 | −1583 | 32/32 |
| 3 | +3.88e7 | +7.98e7 | 2712/2700 |
| 4 | −1.6e11 | −4.8e10 | 3050/3050 |

At depth 3, female has ~2× positive product mass; depth-2 already favors
female by 13× on net P1-source access. Depth-4 products are dominated by
global recurrence and are not interpretable as drive.

---

## 2. Model sensitivity (bounded LIF, female intact)

### 2.1 Targeted relief of the male hop-2 bottleneck

Silence **AN09B017b+c only** (not d, not global params), seeds 11–13:

| Condition | Male P1 | Female P1 | pref |
|---|---|---|---:|
| WT | 0,0,0 | 11,8,8 | −1.00 |
| mAL silence only | 4,5,1 | 18,18,20 | −0.70 |
| bc_sil only | 0,0,1 | 11,8,8 | −0.93 |
| **bc_sil + mAL** | **10,11,3** | **18,18,20** | **−0.42** |
| bcd_sil + mAL (E3) | 9,9,5 | 18,18,20 | −0.42 |

**Male doubles** when the male-specific inhibitory sinks are removed **and**
mAL is off; **female is bit-identical**. This is a local, female-intact
intervention predicted by the hop-2 anatomy — not a global gain scan.

### 2.2 Consistency with earlier grids

- E1: single AN09B017b/c/d silence alone ≈ null (need mAL off to matter).
- E3: combined bcd + mAL was already the best female-intact male drive;
  anatomy now explains **why b and c specifically**.
- E5/E8 amplify both sexes or require structural what-ifs; they do not
  replace this causal pair (mAL + b/c).

---

## 3. Biological hypothesis (not established)

Under the **model’s** sign map and uniform GRN encoding:

> LgLG5 reaches the courtship core with less immediate load on Glu AN09B017b/c
> (and thus less hop-2 inhibition onto P1-related sources), while still
> engaging shared AN09B017/f and downstream excitation. LgLG6 is wired
> preferentially into b/c and IN05B011a, so even when mAL output is blocked,
> residual local inhibition keeps P1 sparse.

This is a **hypothesis about connectome-imposed routing asymmetry**, not a
demonstration that:

- LgLG5/LgLG6 are the physiological pheromone channels in vivo,
- AN09B017b/c are inhibitory onto P1 in the animal (sign is a model assumption
  from Glu→inhibitory default),
- or that relieving b/c produces male-biased mate choice.

Literature (Ryba, Clowney, Kallman) supports mAL as a courtship brake and
vAB3/P1 as a courtship pathway; the **b/c vs g loading asymmetry** is a
connectome-derived prediction that would need receptor-level and functional
tests.

---

## 4. Direct answer

**Where do LgLG6 and LgLG5 first diverge in effective drive to P1?**

**At hop 2.** Both leave the periphery excitatory-only and converge on a
shared AN09B017 hub set. But LgLG6 loads **AN09B017b/c (Glu, P1-inhibitory)**
~170–260× more per cell than LgLG5, producing a 13× more negative net
P1-source access. Female drive is not explained by “less total inhibition at
hop 2” (LgLG5 has plenty of hop-2 inhibition) but by **routing that inhibition
away from the strongest P1-inhibitory types**, leaving AN09B017f and deeper
excitatory paths more effective.

**Local test:** removing b/c transmission with mAL intact does little;
removing b/c **with mAL silenced** raises male P1 ~2× and leaves female
unchanged — the causal signature of a male-specific hop-2 bottleneck.

---

## 5. Claim ledger

| Claim | Kind |
|---|---|
| LgLG6 → AN09B017b/c synapses ≫ LgLG5 | Anatomical |
| Net hop-2 P1-source access 13× more negative for LgLG6 | Anatomical (weighted) |
| bc_sil + mAL raises male P1, female unchanged | Model sensitivity |
| AN09B017b/c inhibit P1 in vivo | **Not established** (sign assumption) |
| This asymmetry causes mate preference | **Not claimed** |
