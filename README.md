# Gay — MaleCNS courtship-circuit orientation package

Bounded LIF + MaleCNS v1.0 connectome experiments on male- vs female-cue
drive to P1/pC1. Preference score `(M−F)/(M+F+ε)` is a **modeled response
index**, not behavior.

**Script map:** [SCRIPTS.md](SCRIPTS.md) · **Architecture synthesis:**
[ARCHITECTURE_SELECTIVITY.md](ARCHITECTURE_SELECTIVITY.md)

---

## Final task goals

1. **Locate** where male sensory drive is gated onto the courtship readout
   (delivered currents onto P1, not source spike counts).
2. **Separate** anatomical fact, model assumption, and biological hypothesis
   — especially AN09B017b/c sign.
3. **Test** female-intact local mechanisms (brake relief, excitatory relay,
   modulatory gates, reward plasticity) without cutting the female path or
   driving P1 directly.
4. **Assess robustness** of any male-release cut-set to transmitter/sign
   uncertainty on b/c.
5. **Synthesize** whether sex selectivity = mAL male-specific I + AN09B017c
   + FLA female-biased E, and what is still unproven in vivo.

**Non-goals:** expanding parameter scans after seeing outcomes; claiming
mate preference; retuning LIF constants; editing `experiment/` upstream.

---

## Known conclusions (in-model)

| # | Claim | Evidence |
|---|---|---|
| 1 | WT male P1 = 0; female = 11,8,8 (seeds 11–13) | `RESULTS.md` |
| 2 | Male signal reaches LgLG6 and AN09B017; dies at P1 E/I veto | `MALE_SIGNAL_AUTOPSY.md` |
| 3 | Last-hop male-specific I: **mAL_m8** (sign-robust), then c/m2b/VES022 | `P1_CURRENT_REPORT.md`, `SIGN_ROBUST_REPORT.md` |
| 4 | Last-hop female-biased E: **FLA001m / FLA003m** | same |
| 5 | AN09B017c loading is male-biased; **inhibitory gate only if c=−1** | `BC_SIGN_ADJUDICATION.md`, `SIGN_ROBUST_REPORT.md` |
| 6 | b/c=+1 ⇒ WT already male>female (phenotype C) | `BC_HOP3_REPORT.md` |
| 7 | No all-sign-map female-intact cut-set without c | `SIGN_ROBUST_REPORT.md` |
| 8 | No hop-1 male∩PAM on courtship core; LgLG→KC = 0 | `CONVERGENCE_REPORT.md` |
| 9 | MB reward learning cannot form male bias | `REWARD_RESULTS.md` |
| 10 | OA / VES022 / SIP106m drive does not release male P1 | `OA_GATE_RESULTS.md` |
| 11 | Architecture: distributed last-hop E/I code + fragile c hinge | `ARCHITECTURE_SELECTIVITY.md` |

---

## Setup

```sh
uv sync
uv run python prepare_data.py .experiment-data   # MaleCNS v1.0, several GB
uv run python run_baselines.py .experiment-data
uv run pytest -q
```

Intervention API: `orientation/` (`output_silence`, `output_gain`,
activation). Seeds 11–13. Do not retune frozen LIF constants in
`orientation/simulate.py`.

---

## Evidence hierarchy

```text
Anatomical (MaleCNS edges, types, NT consensus)
    → Model assumption (Glu→−1, vAB3e/f +1, uniform Poisson)
        → Delivered-current LIF (this package)
            → Behavioral / in vivo  (NOT tested here)
```

Anything about mate choice requires experiments outside this repo.
