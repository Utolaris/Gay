# Gay

Independent Intervention API for male-biased courtship-circuit experiments on
the full MaleCNS v1.0 bounded LIF network (166,606 neurons).

Pre-registered baselines (WT vs mAL output silence) reproduce upstream
follow-up P1 counts. A post-hoc search for conditions that make male P1 >
female P1 is documented in [EXPLORATORY.md](EXPLORATORY.md); it does not
replace those baselines.

Mechanism-search audit (model incompleteness, signed male→P1 pathways,
phenotype A/B/C split, next experiments that keep the female pathway intact)
is in [MECHANISM_SEARCH.md](MECHANISM_SEARCH.md).

E1 brake lesion map (female intact): [E1_PROTOCOL.md](E1_PROTOCOL.md),
[E1_RESULTS.md](E1_RESULTS.md), `e1-results.json`.

E2 male excitatory relay amplification: [E2_PROTOCOL.md](E2_PROTOCOL.md),
[E2_RESULTS.md](E2_RESULTS.md), `e2-results.json`.

## Setup

```sh
uv sync
```

### Prepare the connectome graph

Download public MaleCNS v1.0 tables (~1.1 GB), verify checksums, and build the
classified-neuron graph into `.experiment-data/` (gitignored):

```sh
uv run python prepare_data.py .experiment-data
```

Requires several GB of disk and memory. Source data: MaleCNS v1.0, CC BY 4.0.

## Pre-registered baselines

```sh
uv run python run_baselines.py .experiment-data
# optional: --workers 1|N  (default auto, capped at 4 on this workload)
```

Writes `baseline-results.json`:

- P1/pC1 spikes and Hz
- broad P1-related / mAL / global activity
- matched event hashes
- preference score `(M−F)/(M+F+ε)`

### Trial-level parallelism

`orientation/parallel.py` runs independent (intervention × input × seed) trials
in a process pool. The single-trial Numba kernel and seed semantics are
unchanged. Graph CSR arrays are exported once to
`.experiment-data/graph_mmap/*.npy` and opened `mmap_mode='r'` in workers.

```sh
uv run python regress_parallel.py .experiment-data --workers 4  # must PASS
uv run python benchmark_parallel.py .experiment-data            # 1/2/4/6/8
```

See [PERFORMANCE.md](PERFORMANCE.md).

## Tests

```sh
uv run pytest -q
```

## Intervention API

```python
from orientation import (
    none, output_silence, output_gain, activate_tonic, activate_poisson,
    load_prepared, run_trial, response_metrics, preference_score,
)

net = load_prepared(".experiment-data")
silence = output_silence(net.groups["mAL"])
result = run_trial(net, "candidate_male", seed=11, intervention=silence)
metrics = response_metrics(result, net, silence)
```

Kinds:

- `output_silence` — block outgoing transmission (cells may still spike)
- `output_gain` — multiply outgoing weights by a bounded gain in `[0, max_gain]`
- `activation` — tonic voltage and/or class Poisson drive (not applied to P1 in baselines)

See [PROTOCOL.md](PROTOCOL.md) for the frozen baseline design and
[RESULTS.md](RESULTS.md) for measured outcomes.
