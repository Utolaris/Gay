# Script inventory

Purpose of every remaining script, and how it maps to conclusions.
Obsolete one-off grids (E1–E3_E8, post-hoc preference search, early path
dumps) were removed; git history still has them.

---

## Data & core library

| Path | Role |
|---|---|
| `prepare_data.py` | Download MaleCNS v1.0 tables, build classified-neuron graph into `.experiment-data/` |
| `orientation/simulate.py` | Frozen bounded LIF kernel, Intervention packing, baselines |
| `orientation/intervention.py` | Class-level `output_silence` / `output_gain` / activation |
| `orientation/experimental.py` | Per-neuron gain/tonic/arousal kernel (E-block legacy + learning) |
| `orientation/parallel.py` | Trial-level process pool + graph mmap |
| `orientation/learning.py` | DA-gated KC→MBON plasticity (REWARD block) |
| `orientation/p1_currents.py` | Delivered E/I flux traces onto the 8 P1 cells |
| `run_baselines.py` | Pre-registered WT vs mAL-silence baselines |
| `regress_parallel.py` | Parallel vs serial regression check |
| `tests/test_intervention.py` | Unit tests for intervention API |

---

## Anatomy & pathway (static)

| Script | Output | Conclusion |
|---|---|---|
| `analyze_lglg_pathways.py` | `pathway-decomp-lg6-lg5.json` | Hop-by-hop LgLG6 vs LgLG5; hop-2 male load onto AN09B017b/c |
| `analyze_bc_sign_hop3.py` | `bc-hop3-evidence.json` | b/c annotation/NT evidence; hop-3 roles; sign variants of path mass |
| `analyze_convergence.py` | `convergence-report.json` | No hop-1 male∩PAM on courtship core; OA token on AN09B017; VES022/SIP106m hop-2 candidates |

---

## Mechanism blocks (simulation)

| Script | Protocol / report | Raw | Finding |
|---|---|---|---|
| `run_bc_hop3.py` | `BC_HOP3_*` | `bc-hop3-results.json` | Sign flip b/c=+1 → phenotype C; local bc+mAL raises male, female intact |
| `run_reward.py` | `REWARD_*` | `reward-results.json` | MB reward learning fails (no LgLG→KC); mb_access writes weights, no P1 |
| `run_oa_gate.py` | `OA_GATE_*` | `oa-gate-results.json` | OA / VES022 / SIP106m do not release male P1 |
| `analyze_male_autopsy.py` | `MALE_SIGNAL_AUTOPSY.md` | `male-signal-autopsy.json` | Male signal survives to P1 inputs; cancelled by I>E balance |
| `run_p1_current.py` | `P1_CURRENT_*` | `p1-current-results.json` | Delivered-current ranking; min cut-set mAL_m8+c (map A) |
| `run_sign_robust.py` | `SIGN_ROBUST_*` | `sign-robust-results.json` | mAL_m8 sign-robust; c gate sign-dependent; no all-map cut-set |
| `export_bc_viz.py` | `bc-viz/` | `bc-viz/public/data/trial.json` | Three.js activity export for b/c stimulation |

---

## Documentation map

| Document | Contents |
|---|---|
| `PROTOCOL.md` / `RESULTS.md` | Frozen baseline design and WT/mAL outcomes |
| `MECHANISM_SEARCH.md` | Model incompleteness audit; phenotype A/B/C |
| `LGLG_PATHWAY_REPORT.md` | Why LgLG5 effective drive > LgLG6 at hop-2/3 |
| `BC_SIGN_ADJUDICATION.md` | Is b/c excitatory or inhibitory? Evidence ranking |
| `CONVERGENCE_REPORT.md` | Where male sensory meets modulatory input |
| `MALE_SIGNAL_AUTOPSY.md` | Where male drive is cancelled |
| `P1_CURRENT_REPORT.md` | Last-hop I ranking + minimal cut-set |
| `SIGN_ROBUST_REPORT.md` | Cut-set vs b/c sign maps A/B/C |
| `ARCHITECTURE_SELECTIVITY.md` | Synthesis: mAL + c + FLA as sex-selectivity code |
| `SCRIPTS.md` | This file |

---

## Reproduce pipeline (order)

```sh
uv sync
uv run python prepare_data.py .experiment-data
uv run python run_baselines.py .experiment-data
uv run python analyze_lglg_pathways.py
uv run python analyze_bc_sign_hop3.py
uv run python run_bc_hop3.py .experiment-data --workers 6
uv run python analyze_convergence.py
uv run python run_reward.py .experiment-data --workers 6
uv run python run_oa_gate.py .experiment-data --workers 6
uv run python analyze_male_autopsy.py
uv run python run_p1_current.py .experiment-data
uv run python run_sign_robust.py .experiment-data
uv run pytest -q
```

Optional viewer: `export_bc_viz.py` then `cd bc-viz && bun install && bun run dev`.
