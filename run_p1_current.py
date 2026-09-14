"""Run P1_CURRENT autopsy (see P1_CURRENT_PROTOCOL.md)."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pyarrow.feather as feather

from orientation.p1_currents import build_p1_map, run_p1_current_trial
from orientation.simulate import THRESH_MV, apply_sign_overrides, load_prepared, preference_score

SEEDS = (11, 12, 13)
EPS = 1e-9
MIN_I = 50.0  # pre-registered weight-unit floor for ranking


def type_indices(net, data_dir: Path, types) -> np.ndarray:
    ann = feather.read_table(Path(data_dir) / "annotations.feather").to_pylist()
    by_id = {r["bodyId"]: r for r in ann}
    want = set(types)
    return np.asarray(
        [i for i, bid in enumerate(net.ids) if (by_id.get(int(bid)) or {}).get("type") in want],
        dtype=np.int32,
    )


def rank_sources(trials_m: list[dict], trials_f: list[dict], type_names: list[str], is_inh) -> list[dict]:
    n_types = len(type_names)
    im = np.zeros(n_types)
    iff = np.zeros(n_types)
    em = np.zeros(n_types)
    ef = np.zeros(n_types)
    peak_im = np.zeros(n_types)
    peak_if = np.zeros(n_types)
    # first arrival bin among trials
    first_i_m = np.full(n_types, 10**9)
    first_i_f = np.full(n_types, 10**9)

    for tr in trials_m:
        ib = np.asarray(tr["i_by_type"], dtype=np.float64)  # bins × types
        eb = np.asarray(tr["e_by_type"], dtype=np.float64)
        im += ib.sum(axis=0)
        em += eb.sum(axis=0)
        peak_im = np.maximum(peak_im, ib.max(axis=0))
        for tid in range(n_types):
            col = ib[:, tid]
            nz = np.flatnonzero(col > 0)
            if nz.size:
                first_i_m[tid] = min(first_i_m[tid], int(nz[0]))
    for tr in trials_f:
        ib = np.asarray(tr["i_by_type"], dtype=np.float64)
        eb = np.asarray(tr["e_by_type"], dtype=np.float64)
        iff += ib.sum(axis=0)
        ef += eb.sum(axis=0)
        peak_if = np.maximum(peak_if, ib.max(axis=0))
        for tid in range(n_types):
            col = ib[:, tid]
            nz = np.flatnonzero(col > 0)
            if nz.size:
                first_i_f[tid] = min(first_i_f[tid], int(nz[0]))

    ns = len(trials_m)
    rows = []
    for tid in range(n_types):
        if type_names[tid] == "__none__":
            continue
        Im, If = im[tid] / ns, iff[tid] / ns
        Em, Ef = em[tid] / ns, ef[tid] / ns
        if Im < MIN_I and If < MIN_I and Em < MIN_I and Ef < MIN_I:
            continue
        score = (Im - If) * Im / (Im + If + EPS) if bool(is_inh[tid]) else 0.0
        rows.append(
            {
                "type": type_names[tid],
                "inhibitory": bool(is_inh[tid]),
                "I_male": float(Im),
                "I_female": float(If),
                "I_delta": float(Im - If),
                "E_male": float(Em),
                "E_female": float(Ef),
                "E_delta": float(Em - Ef),
                "peak_I_male": float(peak_im[tid] / 1.0),  # max over trials, not mean
                "peak_I_female": float(peak_if[tid]),
                "first_I_bin_male": None if first_i_m[tid] > 10**8 else int(first_i_m[tid]),
                "first_I_bin_female": None if first_i_f[tid] > 10**8 else int(first_i_f[tid]),
                "male_specific_score": float(score),
            }
        )
    rows.sort(key=lambda r: -r["male_specific_score"])
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data", type=Path)
    ap.add_argument("--out", type=Path, default=Path("p1-current-results.json"))
    args = ap.parse_args()
    started = time.perf_counter()

    net = load_prepared(args.data)
    signs = apply_sign_overrides(net)
    pmap = build_p1_map(net, args.data)
    print(f"P1 n={pmap.p1.size} types={len(pmap.type_names)}", flush=True)

    # Phase A
    phase_a = {"male": [], "female": []}
    for seed in SEEDS:
        for inp, key in (("candidate_male", "male"), ("candidate_female", "female")):
            print(f"PhaseA {key} seed {seed}", flush=True)
            tr = run_p1_current_trial(net, pmap, signs, inp, seed)
            phase_a[key].append(tr)

    ranking = rank_sources(phase_a["male"], phase_a["female"], pmap.type_names, pmap.is_inh_type)
    inh_rank = [r for r in ranking if r["inhibitory"] and r["I_male"] >= MIN_I]
    exc_rank = [r for r in ranking if (not r["inhibitory"]) and r["E_male"] >= MIN_I]
    print("top inhibitory by male-specific score:", flush=True)
    for r in inh_rank[:12]:
        print(
            f"  {r['type'][:20]:20s} Im={r['I_male']:.0f} If={r['I_female']:.0f} "
            f"score={r['male_specific_score']:.1f}",
            flush=True,
        )
    print("top excitatory by E_male:", flush=True)
    for r in sorted(exc_rank, key=lambda x: -x["E_male"])[:8]:
        print(
            f"  {r['type'][:20]:20s} Em={r['E_male']:.0f} Ef={r['E_female']:.0f} d={r['E_delta']:.0f}",
            flush=True,
        )

    # Phase B ablations from frozen ranking
    top_types = [r["type"] for r in inh_rank[:4]]
    mal_types = sorted({t for t in pmap.type_names if t.startswith("mAL_")})
    ladders = []
    for k in (1, 2, 3):
        if len(top_types) >= k:
            ladders.append((f"top{k}", top_types[:k]))
    if top_types:
        combo = list(top_types[:2])
        for t in mal_types:
            if t not in combo:
                combo.append(t)
                break
        ladders.append(("top2_plus_mAL1", combo))

    phase_b = []
    wt_f = [tr["p1_spikes"] for tr in phase_a["female"]]
    wt_m = [tr["p1_spikes"] for tr in phase_a["male"]]
    print(f"WT P1 M={wt_m} F={wt_f}", flush=True)

    min_cut = None
    for name, types in ladders:
        # also silence full mAL class if name includes mAL
        sil_idx = type_indices(net, args.data, types)
        if "mAL1" in name:
            sil_idx = np.unique(np.concatenate([sil_idx, net.groups["mAL"]]))
        gains = np.ones(net.n, dtype=np.float64)
        gains[sil_idx] = 0.0
        ms, fs = [], []
        for seed in SEEDS:
            tm = run_p1_current_trial(
                net, pmap, signs, "candidate_male", seed, gain=gains
            )
            tf = run_p1_current_trial(
                net, pmap, signs, "candidate_female", seed, gain=gains
            )
            ms.append(tm["p1_spikes"])
            fs.append(tf["p1_spikes"])
        row = {
            "condition": name,
            "types": types,
            "male_P1": ms,
            "female_P1": fs,
            "pref": [
                preference_score(a, b) for a, b in zip(ms, fs)
            ],
            "male_gt0_seeds": sum(1 for x in ms if x > 0),
            "female_near_wt": sum(
                1 for i, x in enumerate(fs) if x >= wt_f[i] - 2
            ),
        }
        phase_b.append(row)
        print(json.dumps(row), flush=True)
        if min_cut is None and row["male_gt0_seeds"] >= 2 and row["female_near_wt"] >= 2:
            min_cut = name

    # compact traces: keep one seed male/female full traces; others summary
    def slim_trial(tr):
        return {
            "input": tr["input"],
            "seed": tr["seed"],
            "p1_spikes": tr["p1_spikes"],
            "p1_cell_spikes": tr["p1_cell_spikes"],
            "first_spike_bin": tr["first_spike_bin"],
            "vm_mean": np.mean(tr["vm"], axis=1).tolist(),
            "vm_max": np.max(tr["vm"], axis=1).tolist(),
            "ge_sum": np.sum(tr["ge"], axis=1).tolist(),
            "hi_sum": np.sum(tr["hi"], axis=1).tolist(),
            "e_flux_sum": np.sum(tr["e_flux"], axis=1).tolist(),
            "i_flux_sum": np.sum(tr["i_flux"], axis=1).tolist(),
        }

    report = {
        "kind": "P1_CURRENT_autopsy",
        "protocol_sha256": hashlib.sha256(Path("P1_CURRENT_PROTOCOL.md").read_bytes()).hexdigest(),
        "seeds": list(SEEDS),
        "min_I": MIN_I,
        "WT_male_P1": wt_m,
        "WT_female_P1": wt_f,
        "phase_a_full_seed11": {
            "male": phase_a["male"][0],
            "female": phase_a["female"][0],
        },
        "phase_a_summary": {
            "male": [slim_trial(t) for t in phase_a["male"]],
            "female": [slim_trial(t) for t in phase_a["female"]],
        },
        "ranking_inhibitory": inh_rank[:40],
        "ranking_excitatory": sorted(exc_rank, key=lambda r: -r["E_male"])[:20],
        "phase_b": phase_b,
        "min_cut_set": min_cut,
        "seconds": time.perf_counter() - started,
    }
    args.out.write_text(json.dumps(report))
    print("COMPLETE", report["seconds"], "min_cut", min_cut, flush=True)


if __name__ == "__main__":
    main()
