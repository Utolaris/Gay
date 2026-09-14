# SIGN_ROBUST protocol

Frozen before outcomes. MaleCNS topology, sensory encoding, LIF, P1 IDs,
female pathway unchanged. No P1 drive, no global gain, no deleting all
inhibition, no post-hoc strength tuning, no hand-picked combos that force
AN09B017c.

## Sign maps (only b/c signs differ)

| Map | AN09B017b | AN09B017c | Rationale |
|---|---:|---:|---|
| A | −1 | −1 | package baseline (Glu→−1) |
| B | +1 | +1 | adjudicated / predicted-ACh / vAB3-class |
| C | 0 | 0 | uncertainty control (no fast output) |

All other overrides (vAB3e/f +1, female sensory +1) identical.

## Phase A — per-map delivered currents

For each map × {male, female} × seeds 11,12,13:

- Rebuild signs from scratch (no reuse of another map’s ranking).
- Record per P1 cell: Vm(t), ge, hi, E flux, I flux, first-spike bin.
- Record per presynaptic type: delivered E/I weight flux time course.
- Track named types: mAL_m8, mAL_m2b, AN09B017c, VES022, PVLP048,
  FLA001m, FLA003m.

Rankings (per map, from that map only):

- Male-specific inhibitory score:
  `(I_m − I_f) · I_m / (I_m + I_f + ε)` for types with
  `sign=−1` under **that map** and `I_m ≥ 50`.
- Female-biased excitatory: sort by `E_f − E_m` among `sign=+1` types
  with `E ≥ 50`.

## Phase B — per-map minimal cut-set

Ladder from **that map’s** inhibitory ranking only:

- top1, top2, top3, top4 (types that are inhibitory under that map).
- Ablation: gain=0 (block delivery; source may still spike).
- Primary: male P1 > 0 on ≥2/3 seeds **and** female P1 ≥ map_WT_F − 2
  on ≥2/3 seeds.
- Secondary: male P1 > 0 on 3/3 and female bit-identical to map WT.

Min-cut = first ladder step meeting primary, else best secondary, else null.

## Cross-map robustness labels

| Label | Rule |
|---|---|
| **sign-robust male-expression cut-set** | A set S of types such that for **every** map in {A,B,C}, ablating S (or the map-valid subset S∩inh(map)) meets primary criterion **and** S is built only from types inhibitory in **all** maps (e.g. GABA types) |
| **sign-dependent mechanism** | Works only on some maps (e.g. requires c=−1) |
| **architecture-level** | What holds regardless of cut-set success |

Construction of robust set: intersection of top-k inhibitory types across
maps that remain inhibitory in all maps, then test that intersection as a
cut-set on each map. Do not manually insert c.

## Outputs

- `sign-robust-results.json`
- `SIGN_ROBUST_REPORT.md`
