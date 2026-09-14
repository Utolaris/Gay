# P1_CURRENT protocol

Frozen before ablation outcomes. Topology, sign map, LIF, sensory encoding,
P1 readout IDs, mAL/b/c definitions unchanged. No direct P1 drive. No
global gain. No post-hoc η/threshold retuning. No deleting all inhibition.

## Phase A — instrumented currents (no intervention)

Male and female cue × seeds 11,12,13. Record per P1 cell (8 cells):

- Vm(t) in 10 ms bins
- delivered excitatory weight flux (Σ w of E arrivals in bin)
- delivered inhibitory weight flux (Σ w of I arrivals in bin)
- conductance proxies ge(t), hi(t)
- net flux E−I
- threshold margin = Vm − THRESH_MV
- first-spike latency (ms, or null)

Aggregate every P1 **presynaptic type** (by MaleCNS `type`):

- ∫ delivered I flux dt (sum over bins)
- peak I flux
- first I arrival (ms)
- same for E
- male vs female means over seeds
- male-specific inhibition score:
  `score = (I_m − I_f) · I_m / (I_m + I_f + ε)`
  (large only if type delivers strong I and more to male than female)

Primary ranking key: `score` descending, restricted to `I_m ≥ 50`
(weight units) so tiny edges cannot win.

## Phase B — minimal causal ablation (pre-registered ladder)

From Phase A ranking of **inhibitory** types:

| ID | Ablation (`output_silence` / gain=0 on type) |
|---|---|
| top1 | highest score type |
| top2 | top 1+2 |
| top3 | top 1+2+3 |
| top1+2+mAL | if mAL not in top2, add full mAL class |
| min_cut | smallest k∈{1,2,3,4} with male P1>0 on ≥2/3 seeds **and** female P1 ≥ WT_F − 2 on ≥2/3 seeds |

Ablation blocks **delivered current** only (sources may still spike).
Female path never silenced. All Phase B rows reported.

## Success

- Male P1 > 0 with female near WT ⇒ candidate male-specific last-hop gate.
- No small k works ⇒ sex selectivity is distributed E/I population code.

## Outputs

- `p1-current-results.json` (traces + rankings + ablations)
- `P1_CURRENT_REPORT.md`
