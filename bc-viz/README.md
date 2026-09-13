# bc-viz

Three.js viewer of MaleCNS activity under direct **AN09B017b/c** stimulation.

## Data

`public/data/trial.json` is produced by (from the Gay package root):

```sh
uv run python export_bc_viz.py
```

Conditions (seed 11, 80 Hz Poisson on b/c, 300 ms):

| key | meaning |
|---|---|
| `wt` | no sensory input, no drive |
| `bc_base` | b/c stimulated under Glu→−1 (package default) |
| `bc_exc` | b/c stimulated with b/c signs forced +1 |

Measured: `bc_base` network spikes ≈ 99, P1 = 0; `bc_exc` network spikes ≈ 5.5×10⁴, P1 = 28.

## Run

```sh
bun install
bun run dev      # http://127.0.0.1:5173
bun run build    # static dist/
bun run preview
```

## Interpretation

This is a **model sensitivity instrument**, not a behavioral assay. The wavefront
under `bc_exc` shows how an excitatory map on the male-biased LgLG6→c labeled
line can recruit the central brain; under the inhibitory map the same stimulus
stays local.
