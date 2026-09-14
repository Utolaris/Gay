# OA-gate results

Protocol: [OA_GATE_PROTOCOL.md](OA_GATE_PROTOCOL.md). Raw: `oa-gate-results.json`.
24/24 runs, 0 error. Drive = tonic 12 mV. b/c, mAL, female path unchanged.

| Condition | Male P1 | Female P1 | pref | C? |
|---|---|---|---:|---|
| WT | 0,0,0 | 11,8,8 | −1.00 | no |
| **OA (all OA*)** | **0,0,0** | **11,8,8** | −1.00 | no |
| **VES022** | **0,0,0** | **0,0,1** | −0.33 | no |
| **SIP106m** | **0,0,0** | **11,8,8** | −1.00 | no |

## Answer

**No.** Under frozen signs and intact mAL, octopamine drive does **not**
make the natural male→courtship pathway influence P1.

1. **Global OA activation is inert** on both cues — bit-identical to WT.
   Token OA synapses onto AN09B017 (2–13) are not a usable gate at this
   drive level.
2. **SIP106m activation is inert** on P1 (male and female match WT).
3. **VES022 activation suppresses female** (11,8,8 → 0,0,1) and leaves
   male at 0. VES022 is not a male-biased release valve here; if anything
   it is a shared/female-side brake under this sign map.

## Relation to convergence anatomy

Anatomy said VES022/SIP106m are the best OA×male hop-2 candidates. Function
says activating them (or the OA population) still cannot push male drive
past the innate P1 gate (mAL + hop-2 inhibition + weak male direct path).

Consistent with S-block and REWARD block: **expression of male→P1 is blocked
downstream**, not for lack of a modulatory touch on the ascending line.

## Claim discipline

| Statement | Status |
|---|---|
| OA tonic 12 mV raises male P1 | Does not hold |
| VES022 is a male-biased OA gate onto P1 | Does not hold (female↓, male flat) |
| OA anatomy ⇒ functional release of male courtship drive | Not supported |
| OA cannot matter at any dose / in vivo | Not claimed |

## Reproduce

```sh
uv run python run_oa_gate.py .experiment-data --workers 6 --out oa-gate-results.json
```
