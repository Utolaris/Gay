"""Run REWARD learning block (see REWARD_PROTOCOL.md)."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np

from orientation.learning import (
    ETA,
    N_EPOCHS,
    PAM_HZ,
    W_FLOOR,
    build_plastic_net,
    class_spikes,
    reset_weights,
    run_learn_trial,
    weight_stats,
)
from orientation.simulate import apply_sign_overrides, load_prepared, preference_score

SEEDS = (11, 12, 13)
INPUTS = ("candidate_male", "candidate_female")


def _limit_threads() -> None:
    try:
        from numba import set_num_threads

        set_num_threads(1)
    except Exception:
        pass


_CACHE: dict[str, Any] = {}


def get_pnet(data_dir: str):
    if data_dir not in _CACHE:
        net = load_prepared(data_dir)
        pnet = build_plastic_net(net, Path(data_dir))
        signs = apply_sign_overrides(net)
        _CACHE[data_dir] = (pnet, signs)
    return _CACHE[data_dir]


def measure(pnet, signs, inp: str, seed: int, **kw) -> dict:
    counts, info = run_learn_trial(pnet, signs, inp, seed, **kw)
    p1 = pnet.net.groups["P1_readout"]
    broad = pnet.net.groups["P1_related"]
    return {
        "input": inp,
        "seed": seed,
        "P1_spikes": class_spikes(counts, p1),
        "broad_P1_spikes": class_spikes(counts, broad),
        "KC_spikes": class_spikes(counts, pnet.kc),
        "MBON_spikes": class_spikes(counts, pnet.mbon),
        "PAM_spikes": class_spikes(counts, pnet.pam),
        "network_spikes": int(counts.sum()),
        **info,
    }


def train_pairing(condition: str, epoch: int, seed: int) -> dict:
    """Return flags for the male/female training trials this epoch."""
    # coin for rand_rew uses seed+epoch
    if condition == "male_rew":
        return {
            "male": {"input": "candidate_male", "da_on": True, "reward_hz": PAM_HZ},
            "female": {"input": "candidate_female", "da_on": False, "reward_hz": 0.0},
        }
    if condition == "no_rew":
        return {
            "male": {"input": "candidate_male", "da_on": False, "reward_hz": 0.0},
            "female": {"input": "candidate_female", "da_on": False, "reward_hz": 0.0},
        }
    if condition == "rand_rew":
        coin = np.random.RandomState(seed * 1000 + epoch).rand() < 0.5
        return {
            "male": {
                "input": "candidate_male",
                "da_on": bool(coin),
                "reward_hz": PAM_HZ if coin else 0.0,
            },
            "female": {"input": "candidate_female", "da_on": False, "reward_hz": 0.0},
        }
    if condition == "fem_rew":
        return {
            "male": {"input": "candidate_male", "da_on": False, "reward_hz": 0.0},
            "female": {"input": "candidate_female", "da_on": True, "reward_hz": PAM_HZ},
        }
    if condition == "plast_off":
        return {
            "male": {"input": "candidate_male", "da_on": True, "reward_hz": PAM_HZ},
            "female": {"input": "candidate_female", "da_on": False, "reward_hz": 0.0},
        }
    if condition == "mb_access":
        return {
            "male": {
                "input": "candidate_male",
                "da_on": True,
                "reward_hz": PAM_HZ,
                "kc_access": True,
            },
            "female": {"input": "candidate_female", "da_on": False, "reward_hz": 0.0},
        }
    raise ValueError(condition)


CONDITIONS = (
    "male_rew",
    "no_rew",
    "rand_rew",
    "fem_rew",
    "plast_off",
    "mb_access",
)


def run_one(data_dir: str, condition: str, seed: int) -> dict:
    _limit_threads()
    pnet, signs = get_pnet(data_dir)
    reset_weights(pnet)
    learn_on = condition != "plast_off"

    pre = {inp: measure(pnet, signs, inp, seed, learn_on=False, da_on=False) for inp in INPUTS}
    w_pre = weight_stats(pnet)

    epochs = []
    for ep in range(N_EPOCHS):
        plan = train_pairing(condition, ep, seed)
        rows = {}
        for key in ("male", "female"):
            spec = plan[key]
            rows[key] = measure(
                pnet,
                signs,
                spec["input"],
                seed,
                da_on=bool(spec["da_on"]),
                learn_on=learn_on,
                reward_hz=spec["reward_hz"],
                kc_access=bool(spec.get("kc_access", False)),
            )
        epochs.append(
            {
                "epoch": ep,
                "male": rows["male"],
                "female": rows["female"],
                "weights": weight_stats(pnet),
            }
        )

    post = {inp: measure(pnet, signs, inp, seed, learn_on=False, da_on=False) for inp in INPUTS}
    w_post = weight_stats(pnet)

    # reversal: 20 epochs female+reward regardless (pre-registered)
    rev_epochs = []
    for ep in range(N_EPOCHS):
        rows = {
            "male": measure(
                pnet, signs, "candidate_male", seed, da_on=False, learn_on=learn_on, reward_hz=0.0
            ),
            "female": measure(
                pnet,
                signs,
                "candidate_female",
                seed,
                da_on=learn_on,
                learn_on=learn_on,
                reward_hz=PAM_HZ,
            ),
        }
        rev_epochs.append(
            {"epoch": ep, "male": rows["male"], "female": rows["female"], "weights": weight_stats(pnet)}
        )
    post2 = {inp: measure(pnet, signs, inp, seed, learn_on=False, da_on=False) for inp in INPUTS}
    w_post2 = weight_stats(pnet)

    def pref(block: dict) -> float:
        return preference_score(block["candidate_male"]["P1_spikes"], block["candidate_female"]["P1_spikes"])

    return {
        "condition": condition,
        "seed": seed,
        "eta": ETA,
        "w_floor": W_FLOOR,
        "n_epochs": N_EPOCHS,
        "pre": pre,
        "pre_pref": pref(pre),
        "weights_pre": w_pre,
        "train": epochs,
        "post": post,
        "post_pref": pref(post),
        "weights_post": w_post,
        "reversal_train": rev_epochs,
        "post2": post2,
        "post2_pref": pref(post2),
        "weights_post2": w_post2,
        "delta_P1_male": post["candidate_male"]["P1_spikes"] - pre["candidate_male"]["P1_spikes"],
        "delta_P1_female": post["candidate_female"]["P1_spikes"] - pre["candidate_female"]["P1_spikes"],
        "delta_MbON_male": post["candidate_male"]["MBON_spikes"] - pre["candidate_male"]["MBON_spikes"],
        "delta_KC_male": post["candidate_male"]["KC_spikes"] - pre["candidate_male"]["KC_spikes"],
    }


def summarize(rows: list[dict]) -> list[dict]:
    by: dict[str, list] = {}
    for r in rows:
        by.setdefault(r["condition"], []).append(r)
    out = []
    for cond, items in sorted(by.items()):
        items = sorted(items, key=lambda r: r["seed"])
        out.append(
            {
                "condition": cond,
                "seeds": [r["seed"] for r in items],
                "pre_P1_M": [r["pre"]["candidate_male"]["P1_spikes"] for r in items],
                "pre_P1_F": [r["pre"]["candidate_female"]["P1_spikes"] for r in items],
                "post_P1_M": [r["post"]["candidate_male"]["P1_spikes"] for r in items],
                "post_P1_F": [r["post"]["candidate_female"]["P1_spikes"] for r in items],
                "post2_P1_M": [r["post2"]["candidate_male"]["P1_spikes"] for r in items],
                "post2_P1_F": [r["post2"]["candidate_female"]["P1_spikes"] for r in items],
                "pre_KC_M": [r["pre"]["candidate_male"]["KC_spikes"] for r in items],
                "post_KC_M": [r["post"]["candidate_male"]["KC_spikes"] for r in items],
                "pre_MBON_M": [r["pre"]["candidate_male"]["MBON_spikes"] for r in items],
                "post_MBON_M": [r["post"]["candidate_male"]["MBON_spikes"] for r in items],
                "mean_ratio_post": float(np.mean([r["weights_post"]["mean_ratio"] for r in items])),
                "frac_depressed_post": float(
                    np.mean([r["weights_post"]["frac_depressed"] for r in items])
                ),
                "post_pref": [round(r["post_pref"], 3) for r in items],
                "post2_pref": [round(r["post2_pref"], 3) for r in items],
                "C_post_all_seeds": all(
                    r["post"]["candidate_male"]["P1_spikes"]
                    > r["post"]["candidate_female"]["P1_spikes"]
                    for r in items
                ),
                "weights_moved": any(r["weights_post"]["frac_depressed"] > 0 for r in items),
                "KC_driven_pre": any(r["pre"]["candidate_male"]["KC_spikes"] > 0 for r in items),
            }
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data", type=Path)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out", type=Path, default=Path("reward-results.json"))
    args = ap.parse_args()
    started = time.perf_counter()
    jobs = [(str(args.data), c, s) for c in CONDITIONS for s in SEEDS]
    print(f"jobs={len(jobs)} epochs={N_EPOCHS} eta={ETA}", flush=True)
    rows: list[dict | None] = [None] * len(jobs)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(run_one, *j): i for i, j in enumerate(jobs)}
        done = 0
        for fut in as_completed(futs):
            i = futs[fut]
            try:
                rows[i] = fut.result()
            except Exception as exc:
                rows[i] = {
                    "condition": jobs[i][1],
                    "seed": jobs[i][2],
                    "error": repr(exc),
                }
            done += 1
            print(f"completed {done}/{len(jobs)}", flush=True)
    ok = [r for r in rows if r is not None and "error" not in r]
    summary = summarize(ok) if ok else []
    report = {
        "kind": "REWARD_learning_block",
        "protocol_sha256": hashlib.sha256(Path("REWARD_PROTOCOL.md").read_bytes()).hexdigest(),
        "seeds": list(SEEDS),
        "n_epochs": N_EPOCHS,
        "eta": ETA,
        "w_floor": W_FLOOR,
        "n_jobs": len(jobs),
        "n_ok": len(ok),
        "n_error": len(rows) - len(ok),
        "results": rows,
        "summary": summary,
        "seconds": time.perf_counter() - started,
    }
    args.out.write_text(json.dumps(report, indent=2))
    print("=== SUMMARY ===", flush=True)
    for s in summary:
        print(json.dumps(s), flush=True)
    print("COMPLETE", report["seconds"], flush=True)


if __name__ == "__main__":
    main()
