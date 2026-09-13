# Orientation package conventions

- Frozen LIF constants live in `orientation/simulate.py`. Do not retune them to obtain a target result.
- Interventions are class-level only: output silence, bounded output gain, activation. Never target the P1 readout in sensory baselines.
- Male / female / no-input must share identical network parameters; external schedules are generated once per (input, seed) and reused.
- Seeds for baselines are fixed in advance: 11, 12, 13.
- Preference score is `(M−F)/(M+F+eps)` on P1 spikes. It is a modeled response index, not behavior.
- Report all runs, including null or female-biased outcomes.
- Upstream `experiment/` and `experiment/followup/` are read-only references in this package's workflow; do not edit them here.
- Python tooling: `uv` in this directory. Frontend remains Bun at the repo root and is out of scope for orientation.
