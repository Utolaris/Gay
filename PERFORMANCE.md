# Performance: trial-level parallelism and prange evaluation

Model outputs are unchanged. Production `_simulate` is not modified.

## Trial-level multiprocessing

| Component | Design |
|---|---|
| Parallel axis | independent trials (intervention × input × seed) |
| Default workers | `min(CPU, n_trials, 4)` — auto; override `--workers` |
| Kernel / seeds | identical to serial (`run_trial` + `sensory_event_schedule`) |
| Graph sharing | once-written `graph_mmap/*.npy`, `mmap_mode='r'` per worker |
| Numba threads | forced to 1 inside workers (avoid oversubscription) |

### Regression (required)

18 baseline trials, serial vs 4 workers:

```text
match: true
mismatches: []
serial 56.4 s → parallel 24.0 s (2.35×)
```

Reproduce: `uv run python regress_parallel.py .experiment-data --workers 4`

### Benchmark (8-core Apple silicon, 18 trials)

| Workers | Wall (s) | CPU (s) | CPU util | Child peak RSS (MB) | Speedup |
|---:|---:|---:|---:|---:|---:|
| 1 | 55.2 | 55.2 | 1.00 | 1.5 | 1.00× |
| 2 | 31.9 | 65.4 | 2.05 | 317 | 1.73× |
| **4** | **25.3** | **94.5** | **3.73** | **319** | **2.18×** |
| 6 | 33.1 | 186.6 | 5.63 | 319 | 1.67× |
| 8 | 38.5 | 203.6 | 5.29 | 319 | 1.44× |

**Best wall time at 4 workers** for an 18-trial batch. 6–8 oversubscribe
(pool startup + uneven last-wave). Child RSS ~319 MB indicates shared mmap
pages rather than a full private graph copy per worker.

Larger grids (E1 90 trials, E2 60 trials) benefit more from higher worker
counts; pass `--workers` explicitly when `n_trials` is large.

## `_simulate` `prange` evaluation (not adopted)

Probe: `eval_prange.py` — parallelizes **only** the per-neuron subthreshold
membrane loop (`prange`); synaptic scatter stays **serial** (shared-write race
on `pending_e`/`pending_i`).

| Numba threads | Spike counts vs prod | Extrema vs prod | Speedup |
|---:|---|---|---:|
| 1 | equal | equal | 1.02× |
| 2 | equal | differ* | 1.83× |
| 4 | equal | differ* | 2.84× |
| 8 | equal | differ* | 2.89× |

\* Probe tracks `vmax` after external drive (production tracks after membrane
update). Spike trains are identical; extrema are diagnostics only.

**Decision:** do **not** change production `_simulate` in this pass.

- Trial-level process parallelism already uses multiple cores without touching
  dynamics.
- Combining process workers × `prange` threads would oversubscribe.
- A future single-process mode could adopt membrane-only `prange` if a full
  multi-seed regression shows identical counts **and** matched extrema scan.

Reproduce: `uv run python eval_prange.py .experiment-data`
