# Where the male signal dies

Static LIF autopsy, male cue, seed 11. Script: `analyze_male_autopsy.py`.
Raw: `male-signal-autopsy.json`. Signs, mAL, b/c, female path unchanged.

## Answer in one line

**The male signal is not eliminated on the ascending line.** LgLG6 and
AN09B017 fire. It is **cancelled at the last hop onto P1**, where ~470
inhibitory sources outweigh ~627 excitatory sources. mAL is the largest
named single brake, but **not the only one**.

---

## Spike table (response window, male cue)

| Condition | LgLG6 | b | c | e | f | AN05B035 | AN03A008 | mAL | P1_exc_src | P1_inh_src | **P1** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| WT | 537 | 2 | 36 | 21 | 50 | 82 | 15 | 77 | **255** | **296** | **0** |
| mAL_sil | 534 | 2 | 33 | 28 | 50 | 82 | 15 | 115* | 314 | 382 | 4 |
| bc_sil | 543 | 3 | 35 | 24 | 49 | 83 | 16 | 82 | 286 | 329 | 0 |
| **bc+mAL** | 533 | 3 | 34 | 29 | 49 | 84 | 14 | 119* | 388 | 465 | **10** |
| mAL+IN11a | 606 | 3 | 47 | 26 | 53 | 85 | 19 | 112* | 324 | 403 | 1 |
| mAL+VES022 | 532 | 2 | 34 | 28 | 50 | 81 | 15 | 121* | 320 | 395 | 4 |
| mAL+a008 | 534 | 2 | 36 | 28 | 49 | 80 | 36 | 128* | 378 | 470 | 1 |
| **mAL+all P1-inh** | 532 | 3 | 36 | 28 | 52 | 83 | 24 | — | **7640** | 7478 | **267** |

\*mAL still *spikes* under output_silence (transmission blocked, not cell silenced).

Female P1 under the same grid stays ~11–19 unless P1-inh is globally cut.

---

## Path autopsy

```text
LgLG6 (537 spikes)          ← signal PRESENT
    ↓
AN09B017b/c/e/f, AN05B035   ← signal PRESENT (hub fires)
    ↓
P1_exc_src 255  vs  P1_inh_src 296   ← BALANCE TIPS TO INHIBITION
    ↓
P1 = 0                      ← signal CANCELLED HERE
```

| Stage | Status |
|---|---|
| Sensory GRN | Alive (LgLG6 537) |
| AN09B017 hub | Alive (c=36, f=50, e=21) |
| Shared relays AN05B035 / AN03A008 | Alive |
| **P1 input E/I balance** | **Inhibition wins (296>255)** |
| P1 output | 0 |

So “消灭” is **not** a single labeled-line cutoff (not LgLG, not AN09B017
silence). It is **last-stage inhibitory veto** on the P1 microcircuit.

---

## Which brakes matter

| Removal | ΔP1 (male, seed 11) | Notes |
|---|---:|---|
| mAL only | 0 → **4** | largest named brake; residual inh remains |
| b/c only | 0 → 0 | need mAL off to matter |
| mAL + b/c | 0 → **10** | best local pair (matches S2) |
| mAL + IN11a | 0 → 1 | not a male gate |
| mAL + VES022 | 0 → 4 | no extra male release |
| AN03A008 +12 mV under mAL | 0 → 1 | amplifying male E relay insufficient |
| mAL + **all 470 P1-inh sources** | 0 → **267** | latent excitatory drive is huge |

Interpretation: male drive is **already at the P1 input layer** (exc sources
fire). It fails because **distributed inhibition** (mAL + SIP/VES/LH/… the
470-source set) holds the E/I balance below P1 threshold. Local pair
(mAL+b/c) moves the needle; only wiping essentially all P1 inhibition
unmasks the full latent drive (267 spikes, female 266 — not male-specific).

---

## Relation to earlier blocks

| Block | Consistent? |
|---|---|
| S2: bc+mAL best local male drive | Yes (P1=10 here) |
| E1: no single non-mAL brake | Yes |
| REWARD / OA gate nulls | Yes — they do not rebalance P1 E/I |
| Sign flip bc=+1 → phenotype C | Different lever: turns a hop-2 sink into drive, not a last-hop brake removal |

---

## Claim discipline

| Statement | Status |
|---|---|
| Male cue evokes spikes in LgLG6 and AN09B017 | Holds (model) |
| Male cue fails at P1 because inh sources outfire exc sources | Holds (seed 11) |
| mAL is the unique brake | Does not hold |
| Removing all P1 inhibition is male-specific | Does not hold (F≈266 too) |
| This is the in vivo site of suppression | Model statement only |

## Reproduce

```sh
uv run python analyze_male_autopsy.py
```
