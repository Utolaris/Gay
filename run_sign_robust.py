"""SIGN_ROBUST_CUTSET: male-specific P1 gate vs b/c sign uncertainty."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pyarrow.feather as feather

from orientation.p1_currents import build_p1_map, run_p1_current_trial
from orientation.simulate import (
    apply_sign_overrides,
    load_prepared,
    preference_score,
)

SEEDS = (11, 12, 13)
EPS = 1e-9
MIN_FLUX = 50.0
TRACKED = (
    "mAL_m8",
    "mAL_m2b",
    "AN09B017c",
    "VES022",
    "PVLP048",
    "FLA001m",
    "FLA003m",
)
MAPS = {
    "A_base": {"AN09B017b": -1, "AN09B017c": -1},
    "B_exc": {"AN09B017b": 1, "AN09B017c": 1},
    "C_zero": {"AN09B017b": 0, "AN09B017c": 0},
}


def type_indices(net, data_dir: Path, types) -> np.ndarray:
    ann = feather.read_table(Path(data_dir) / "annotations.feather").to_pylist()
    by_id = {r["bodyId"]: r for r in ann}
    want = set(types)
    return np.asarray(
        [i for i, bid in enumerate(net.ids) if (by_id.get(int(bid)) or {}).get("type") in want],
        dtype=np.int32,
    )


def make_signs(net, base: np.ndarray, overrides: dict[str, int], data_dir: Path) -> np.ndarray:
    signs = base.copy()
    for t, s in overrides.items():
        signs[type_indices(net, data_dir, [t])] = int(s)
    return signs


def type_sign_map(signs: np.ndarray, type_ids: np.ndarray, n_types: int) -> np.ndarray:
    """Sign per type id: first member's sign (all same type should share)."""
    sg = np.zeros(n_types, dtype=np.int8)
    seen = np.zeros(n_types, dtype=np.bool_)
    for i in range(len(type_ids)):
        tid = int(type_ids[i])
        if not seen[tid]:
            seen[tid] = True
            sg[tid] = signs[i]
    return sg


def rank_map(
    trials_m: list[dict],
    trials_f: list[dict],
    type_names: list[str],
    type_signs: np.ndarray,
) -> dict:
    n_types = len(type_names)
    im = np.zeros(n_types)
    iff = np.zeros(n_types)
    em = np.zeros(n_types)
    ef = np.zeros(n_types)
    peak_im = np.zeros(n_types)
    peak_if = np.zeros(n_types)
    first_m = np.full(n_types, 10**9)
    first_f = np.full(n_types, 10**9)

    for tr in trials_m:
        ib = np.asarray(tr["i_by_type"], dtype=np.float64)
        eb = np.asarray(tr["e_by_type"], dtype=np.float64)
        im += ib.sum(axis=0)
        em += eb.sum(axis=0)
        peak_im = np.maximum(peak_im, ib.max(axis=0))
        for tid in range(n_types):
            nz = np.flatnonzero(ib[:, tid] > 0)
            if nz.size:
                first_m[tid] = min(first_m[tid], int(nz[0]))
    for tr in trials_f:
        ib = np.asarray(tr["i_by_type"], dtype=np.float64)
        eb = np.asarray(tr["e_by_type"], dtype=np.float64)
        iff += ib.sum(axis=0)
        ef += eb.sum(axis=0)
        peak_if = np.maximum(peak_if, ib.max(axis=0))
        for tid in range(n_types):
            nz = np.flatnonzero(ib[:, tid] > 0)
            if nz.size:
                first_f[tid] = min(first_f[tid], int(nz[0]))

    ns = len(trials_m)
    inh_rows = []
    exc_rows = []
    for tid in range(n_types):
        name = type_names[tid]
        if name == "__none__":
            continue
        Im, If = im[tid] / ns, iff[tid] / ns
        Em, Ef = em[tid] / ns, ef[tid] / ns
        sg = int(type_signs[tid])
        if sg < 0 and Im >= MIN_FLUX:
            score = (Im - If) * Im / (Im + If + EPS)
            inh_rows.append(
                {
                    "type": name,
                    "I_male": float(Im),
                    "I_female": float(If),
                    "I_delta": float(Im - If),
                    "peak_I_male": float(peak_im[tid]),
                    "peak_I_female": float(peak_if[tid]),
                    "first_I_bin_male": None if first_m[tid] > 10**8 else int(first_m[tid]),
                    "first_I_bin_female": None if first_f[tid] > 10**8 else int(first_f[tid]),
                    "male_specific_score": float(score),
                    "sign": sg,
                }
            )
        if sg > 0 and max(Em, Ef) >= MIN_FLUX:
            exc_rows.append(
                {
                    "type": name,
                    "E_male": float(Em),
                    "E_female": float(Ef),
                    "E_delta": float(Em - Ef),
                    "sign": sg,
                }
            )
    inh_rows.sort(key=lambda r: -r["male_specific_score"])
    exc_rows.sort(key=lambda r: -r["E_delta"])  # most female-biased first
    return {"inhibitory": inh_rows, "excitatory_female_bias": exc_rows}


def tracked_summary(
    trials_m: list[dict],
    trials_f: list[dict],
    type_names: list[str],
    type_signs: np.ndarray,
) -> dict:
    name_to_tid = {n: i for i, n in enumerate(type_names)}
    out = {}
    for name in TRACKED:
        tid = name_to_tid.get(name)
        if tid is None:
            out[name] = None
            continue
        def acc(trials, key):
            tot = peak = 0.0
            first = None
            for tr in trials:
                arr = np.asarray(tr[key], dtype=np.float64)[:, tid]
                tot += float(arr.sum())
                peak = max(peak, float(arr.max()))
                nz = np.flatnonzero(arr > 0)
                if nz.size:
                    b = int(nz[0])
                    first = b if first is None else min(first, b)
            n = max(len(trials), 1)
            return {"mean_total": tot / n, "peak": peak, "first_bin": first}

        out[name] = {
            "sign": int(type_signs[tid]),
            "I_male": acc(trials_m, "i_by_type"),
            "I_female": acc(trials_f, "i_by_type"),
            "E_male": acc(trials_m, "e_by_type"),
            "E_female": acc(trials_f, "e_by_type"),
        }
    return out


def slim_trace(tr: dict) -> dict:
    vm = np.asarray(tr["vm"])
    ge = np.asarray(tr["ge"])
    hi = np.asarray(tr["hi"])
    ef = np.asarray(tr["e_flux"])
    iff = np.asarray(tr["i_flux"])
    netc = ef - iff
    return {
        "input": tr["input"],
        "seed": tr["seed"],
        "p1_spikes": tr["p1_spikes"],
        "p1_cell_spikes": tr["p1_cell_spikes"],
        "first_spike_bin": tr["first_spike_bin"],
        "vm_per_cell_mean_over_time": vm.mean(axis=0).tolist(),
        "vm_max_over_time": vm.max(axis=0).tolist(),
        "threshold_margin_min": (vm - (-45.0)).min(axis=0).tolist(),
        "E_flux_total": ef.sum(axis=0).tolist(),
        "I_flux_total": iff.sum(axis=0).tolist(),
        "net_flux_total": netc.sum(axis=0).tolist(),
        "ge_peak": ge.max(axis=0).tolist(),
        "hi_peak": hi.max(axis=0).tolist(),
    }


def search_cutset(net, pmap, signs, data_dir: Path, inh_rank: list[dict], wt_f: list[int]) -> list[dict]:
    top = [r["type"] for r in inh_rank[:4]]
    rows = []
    for k in (1, 2, 3, 4):
        if len(top) < k:
            break
        types = top[:k]
        idx = type_indices(net, data_dir, types)
        gain = np.ones(net.n, dtype=np.float64)
        gain[idx] = 0.0
        ms, fs = [], []
        for seed in SEEDS:
            tm = run_p1_current_trial(net, pmap, signs, "candidate_male", seed, gain=gain)
            tf = run_p1_current_trial(net, pmap, signs, "candidate_female", seed, gain=gain)
            ms.append(tm["p1_spikes"])
            fs.append(tf["p1_spikes"])
        male_gt0 = sum(1 for x in ms if x > 0)
        female_near = sum(1 for i, x in enumerate(fs) if x >= wt_f[i] - 2)
        female_ident = sum(1 for i, x in enumerate(fs) if x == wt_f[i])
        primary = male_gt0 >= 2 and female_near >= 2
        secondary = male_gt0 == 3 and female_ident == 3
        rows.append(
            {
                "k": k,
                "types": types,
                "male_P1": ms,
                "female_P1": fs,
                "pref": [preference_score(a, b) for a, b in zip(ms, fs)],
                "male_gt0_seeds": male_gt0,
                "female_near_wt": female_near,
                "female_identical": female_ident,
                "primary": primary,
                "secondary": secondary,
            }
        )
        print(f"  cut top{k} {types} M={ms} F={fs} P={primary} S={secondary}", flush=True)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data", type=Path)
    ap.add_argument("--out", type=Path, default=Path("sign-robust-results.json"))
    args = ap.parse_args()
    started = time.perf_counter()

    net = load_prepared(args.data)
    base_signs = apply_sign_overrides(net)
    pmap = build_p1_map(net, args.data)
    report: dict = {
        "kind": "SIGN_ROBUST_CUTSET",
        "protocol_sha256": hashlib.sha256(Path("SIGN_ROBUST_PROTOCOL.md").read_bytes()).hexdigest(),
        "seeds": list(SEEDS),
        "min_flux": MIN_FLUX,
        "maps": {},
    }

    per_map_inh_names: dict[str, list[str]] = {}

    for map_id, overrides in MAPS.items():
        print(f"\n===== MAP {map_id} {overrides} =====", flush=True)
        signs = make_signs(net, base_signs, overrides, args.data)
        tsg = type_sign_map(signs, pmap.type_ids, len(pmap.type_names))

        phase = {"male": [], "female": []}
        for seed in SEEDS:
            for inp, key in (("candidate_male", "male"), ("candidate_female", "female")):
                print(f"  PhaseA {key} seed {seed}", flush=True)
                phase[key].append(run_p1_current_trial(net, pmap, signs, inp, seed))

        wt_m = [t["p1_spikes"] for t in phase["male"]]
        wt_f = [t["p1_spikes"] for t in phase["female"]]
        print(f"  WT P1 M={wt_m} F={wt_f}", flush=True)

        ranks = rank_map(phase["male"], phase["female"], pmap.type_names, tsg)
        print("  top inhibitory:", flush=True)
        for r in ranks["inhibitory"][:10]:
            print(
                f"    {r['type'][:18]:18s} Im={r['I_male']:.0f} If={r['I_female']:.0f} "
                f"score={r['male_specific_score']:.1f} sg={r['sign']}",
                flush=True,
            )
        print("  female-biased E:", flush=True)
        for r in ranks["excitatory_female_bias"][:6]:
            print(
                f"    {r['type'][:18]:18s} Em={r['E_male']:.0f} Ef={r['E_female']:.0f} d={r['E_delta']:.0f}",
                flush=True,
            )

        track = tracked_summary(phase["male"], phase["female"], pmap.type_names, tsg)
        print("  tracked:", flush=True)
        for name, info in track.items():
            if not info:
                print(f"    {name}: ABSENT", flush=True)
                continue
            print(
                f"    {name:12s} sg={info['sign']:+d} "
                f"Im={info['I_male']['mean_total']:.0f} If={info['I_female']['mean_total']:.0f} "
                f"Em={info['E_male']['mean_total']:.0f} Ef={info['E_female']['mean_total']:.0f}",
                flush=True,
            )

        print("  cut-set search:", flush=True)
        cuts = search_cutset(net, pmap, signs, args.data, ranks["inhibitory"], wt_f)
        min_cut = None
        for row in cuts:
            if row["primary"] and min_cut is None:
                min_cut = row
        if min_cut is None:
            for row in cuts:
                if row["secondary"]:
                    min_cut = row
                    break

        per_map_inh_names[map_id] = [r["type"] for r in ranks["inhibitory"][:8]]
        report["maps"][map_id] = {
            "overrides": overrides,
            "WT_male_P1": wt_m,
            "WT_female_P1": wt_f,
            "traces": {
                "male": [slim_trace(t) for t in phase["male"]],
                "female": [slim_trace(t) for t in phase["female"]],
            },
            "full_seed11": {"male": phase["male"][0], "female": phase["female"][0]},
            "ranking_inhibitory": ranks["inhibitory"][:30],
            "ranking_excitatory_female_bias": ranks["excitatory_female_bias"][:15],
            "tracked": track,
            "cutset_ladder": cuts,
            "min_cut": min_cut,
        }

    # Cross-map robust set: types inhibitory in ALL maps and in each map's top list
    sets = []
    for mid in MAPS:
        inh_types = {r["type"] for r in report["maps"][mid]["ranking_inhibitory"]}
        sets.append(inh_types)
    always_inh_top = set.intersection(*sets) if sets else set()
    # rank by mean score across maps among always-inh
    mean_score = {}
    for t in always_inh_top:
        sc = []
        for mid in MAPS:
            for r in report["maps"][mid]["ranking_inhibitory"]:
                if r["type"] == t:
                    sc.append(r["male_specific_score"])
                    break
        mean_score[t] = float(np.mean(sc)) if sc else 0.0
    robust_rank = sorted(mean_score.items(), key=lambda kv: -kv[1])
    report["cross_map"] = {
        "inhibitory_in_all_maps_top30": sorted(always_inh_top),
        "mean_score": mean_score,
        "robust_rank": [{"type": t, "mean_score": s} for t, s in robust_rank],
    }

    # Test intersection top-k as robust cut-set on each map
    robust_candidates = [t for t, _ in robust_rank if t in always_inh_top][:4]
    report["cross_map"]["robust_candidates"] = robust_candidates
    robust_tests = []
    for k in (1, 2, 3, 4):
        if len(robust_candidates) < k:
            break
        types = robust_candidates[:k]
        per_map = {}
        all_primary = True
        all_secondary = True
        for map_id, overrides in MAPS.items():
            signs = make_signs(net, base_signs, overrides, args.data)
            idx = type_indices(net, args.data, types)
            gain = np.ones(net.n, dtype=np.float64)
            gain[idx] = 0.0
            wt_f = report["maps"][map_id]["WT_female_P1"]
            ms, fs = [], []
            for seed in SEEDS:
                tm = run_p1_current_trial(net, pmap, signs, "candidate_male", seed, gain=gain)
                tf = run_p1_current_trial(net, pmap, signs, "candidate_female", seed, gain=gain)
                ms.append(tm["p1_spikes"])
                fs.append(tf["p1_spikes"])
            male_gt0 = sum(1 for x in ms if x > 0)
            female_near = sum(1 for i, x in enumerate(fs) if x >= wt_f[i] - 2)
            female_ident = sum(1 for i, x in enumerate(fs) if x == wt_f[i])
            primary = male_gt0 >= 2 and female_near >= 2
            secondary = male_gt0 == 3 and female_ident == 3
            all_primary = all_primary and primary
            all_secondary = all_secondary and secondary
            per_map[map_id] = {
                "male_P1": ms,
                "female_P1": fs,
                "primary": primary,
                "secondary": secondary,
            }
            print(f"robust top{k} {types} {map_id} M={ms} F={fs} P={primary}", flush=True)
        robust_tests.append(
            {
                "k": k,
                "types": types,
                "per_map": per_map,
                "all_maps_primary": all_primary,
                "all_maps_secondary": all_secondary,
            }
        )
    report["cross_map"]["robust_cutset_tests"] = robust_tests
    label = None
    for row in robust_tests:
        if row["all_maps_primary"]:
            label = {
                "kind": "sign-robust male-expression cut-set" if row["all_maps_secondary"] else "sign-robust primary-only",
                "types": row["types"],
                "k": row["k"],
            }
            break
    report["cross_map"]["label"] = label

    args.out.write_text(json.dumps(report))
    print("\nLABEL", label, flush=True)
    print("COMPLETE", time.perf_counter() - started, flush=True)


if __name__ == "__main__":
    main()
