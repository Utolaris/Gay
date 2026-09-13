# Gay

Independent Intervention API for male-biased courtship-circuit experiments on
the full MaleCNS v1.0 bounded LIF network (166,606 neurons).

Pre-registered baselines (WT vs mAL output silence) reproduce upstream
follow-up P1 counts. A post-hoc search for conditions that make male P1 >
female P1 is documented in [EXPLORATORY.md](EXPLORATORY.md); it does not
replace those baselines.

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
```

Writes `baseline-results.json`:

- P1/pC1 spikes and Hz
- broad P1-related / mAL / global activity
- matched event hashes
- preference score `(M−F)/(M+F+ε)`

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
