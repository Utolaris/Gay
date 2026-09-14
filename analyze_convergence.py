"""Locate natural convergence of male sensory + reward/modulatory onto courtship paths.

Static anatomical analysis. No simulation. Outputs convergence-report.json + printout.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow.feather as feather
from scipy.sparse import load_npz

ROOT = Path(".experiment-data")
P1_IDS = [12442, 16719, 17867, 20117, 20803, 23968, 519518, 522419]


def main() -> None:
    meta = json.loads((ROOT / "full-meta.json").read_text())
    nodes = np.load(ROOT / "full-nodes.npz")
    ids = nodes["ids"]
    ann = feather.read_table(ROOT / "annotations.feather").to_pylist()
    by_id = {r["bodyId"]: r for r in ann}
    W = load_npz(ROOT / "full-graph.npz").tocsr()
    indptr, indices, data = W.indptr, W.indices, W.data
    Wcsc = W.tocsc()

    def ntype(i: int) -> str:
        return (by_id.get(int(ids[i]), {}) or {}).get("type") or "?"

    def outs(i: int):
        s, e = indptr[i], indptr[i + 1]
        return indices[s:e], data[s:e]

    def ins(i: int):
        s, e = Wcsc.indptr[i], Wcsc.indptr[i + 1]
        return Wcsc.indices[s:e], Wcsc.data[s:e]

    p1 = np.searchsorted(ids, P1_IDS)
    p1set = set(int(x) for x in p1)
    male = np.asarray(meta["groups"]["candidate_male"], np.int32)
    female = np.asarray(meta["groups"]["candidate_female"], np.int32)
    mal = np.asarray(meta["groups"]["mAL"], np.int32)
    p1_related = np.asarray(meta["groups"]["P1_related"], np.int32)

    def cells_prefix(prefixes: tuple[str, ...]) -> np.ndarray:
        return np.asarray(
            [
                i
                for i in range(len(ids))
                if any((ntype(i) or "").startswith(p) for p in prefixes)
            ],
            dtype=np.int32,
        )

    pam = cells_prefix(("PAM",))
    ppl = cells_prefix(("PPL",))
    oa = cells_prefix(("OA-", "OAN"))
    da = np.unique(np.concatenate([pam, ppl]))
    # broader modulatory: DA + OA
    mod = np.unique(np.concatenate([da, oa]))
    # courtship core: pC1* + mAL + P1 readout
    pc1 = cells_prefix(("pC1",))
    court = np.unique(np.concatenate([pc1, mal, p1]))

    print("=== populations ===")
    for name, arr in [
        ("male_LgLG", male),
        ("female_LgLG", female),
        ("PAM", pam),
        ("PPL", ppl),
        ("OA", oa),
        ("mod_all", mod),
        ("pC1", pc1),
        ("mAL", mal),
        ("P1_readout", p1),
        ("courtship_core", court),
    ]:
        print(f"{name:16s} n={len(arr)}")

    # accumulate hop-1 weighted maps
    court_set = set(int(x) for x in court)
    male_to = defaultdict(float)
    female_to = defaultdict(float)
    mod_to = defaultdict(float)
    for i in male:
        tg, w = outs(int(i))
        for t, ww in zip(tg, w):
            male_to[int(t)] += float(ww)
    for i in female:
        tg, w = outs(int(i))
        for t, ww in zip(tg, w):
            female_to[int(t)] += float(ww)
    for i in mod:
        tg, w = outs(int(i))
        for t, ww in zip(tg, w):
            mod_to[int(t)] += float(ww)

    print("\n=== hop-1 overlap: male targets ∩ modulatory targets ===")
    male_targets = set(male_to)
    mod_targets = set(mod_to)
    both = male_targets & mod_targets
    print(f"male hop1 targets {len(male_targets)}, mod hop1 targets {len(mod_targets)}, overlap {len(both)}")

    # Among overlap, who also reaches courtship within 2 hops (direct →court or →X→court)
    def reaches_courtship(t: int, max_hops: int = 2) -> tuple[int, float]:
        """Return (min_hops, signed-free product mass of t→…→P1 or court via direct P1)."""
        # direct to P1
        tg, w = outs(t)
        mass = 0.0
        for tt, ww in zip(tg, w):
            if int(tt) in p1set:
                mass += float(ww)
        if mass > 0:
            return 1, mass
        if max_hops < 2:
            return 0, 0.0
        # hop2 via any
        best = 0.0
        for tt, ww in zip(tg, w):
            tg2, w2 = outs(int(tt))
            for t3, ww2 in zip(tg2, w2):
                if int(t3) in p1set:
                    best += float(ww) * float(ww2)
        if best > 0:
            return 2, best
        return 0, 0.0

    # Also score courtship-core membership and output to court
    rows = []
    for t in both:
        mw = male_to[t]
        fw = female_to.get(t, 0.0)
        ow = mod_to[t]
        # breakdown modulatory sources
        pam_w = oa_w = ppl_w = 0.0
        # incoming from mod classes
        pres, pw = ins(t)
        for p, ww in zip(pres, pw):
            p = int(p)
            tp = ntype(p)
            if tp.startswith("PAM"):
                pam_w += float(ww)
            elif tp.startswith("PPL"):
                ppl_w += float(ww)
            elif tp.startswith("OA-") or tp.startswith("OAN"):
                oa_w += float(ww)
        # output to courtship core
        to_court = 0.0
        to_p1 = 0.0
        tg, w = outs(t)
        for tt, ww in zip(tg, w):
            tt = int(tt)
            if tt in p1set:
                to_p1 += float(ww)
            if tt in court_set:
                to_court += float(ww)
        hops, p1mass = reaches_courtship(t, 2)
        rows.append(
            {
                "idx": int(t),
                "id": int(ids[t]),
                "type": ntype(t),
                "from_male": mw,
                "from_female": fw,
                "from_mod": ow,
                "from_PAM": pam_w,
                "from_PPL": ppl_w,
                "from_OA": oa_w,
                "to_P1": to_p1,
                "to_court": to_court,
                "p1_path_hops": hops,
                "p1_path_mass": p1mass,
                "in_court_core": t in court_set,
                "male_bias_in": mw / (fw + 1.0),
            }
        )

    # rank by male*mod and courtship access
    def score(r):
        return r["from_male"] * (r["from_PAM"] + r["from_PPL"] + r["from_OA"] + 1e-9) * (
            r["to_P1"] + 0.01 * r["to_court"] + (r["p1_path_mass"] if r["p1_path_hops"] else 0.0)
        )

    rows.sort(key=score, reverse=True)
    print("\n=== top 25 convergence candidates (male ∩ mod at hop-1, ranked by male×mod×courtship) ===")
    print(
        f"{'type':16s} {'m→T':>7s} {'f→T':>7s} {'PAM':>7s} {'OA':>7s} {'PPL':>6s} "
        f"{'→P1':>6s} {'→court':>7s} {'hops':>4s} {'p1mass':>10s}"
    )
    for r in rows[:25]:
        print(
            f"{r['type'][:16]:16s} {r['from_male']:7.0f} {r['from_female']:7.0f} "
            f"{r['from_PAM']:7.0f} {r['from_OA']:7.0f} {r['from_PPL']:6.0f} "
            f"{r['to_P1']:6.0f} {r['to_court']:7.0f} {r['p1_path_hops']:4d} {r['p1_path_mass']:10.0f}"
        )

    # Restrict to those with courtship access
    court_acc = [r for r in rows if r["to_P1"] > 0 or r["to_court"] > 0 or r["in_court_core"] or r["p1_path_hops"] > 0]
    print(f"\noverlap with courtship access: {len(court_acc)} / {len(rows)}")
    print("\n=== top 20 with courtship access ===")
    for r in court_acc[:20]:
        print(
            f"{r['type'][:18]:18s} id={r['id']} male={r['from_male']:.0f} fem={r['from_female']:.0f} "
            f"PAM={r['from_PAM']:.0f} OA={r['from_OA']:.0f} →P1={r['to_P1']:.0f} →court={r['to_court']:.0f} "
            f"core={r['in_court_core']}"
        )

    # Type-level aggregation for courtship-access cells
    by_type = defaultdict(lambda: defaultdict(float))
    type_n = defaultdict(int)
    for r in court_acc:
        t = r["type"]
        type_n[t] += 1
        for k in ("from_male", "from_female", "from_PAM", "from_OA", "from_PPL", "to_P1", "to_court", "p1_path_mass"):
            by_type[t][k] += r[k]
    print("\n=== type-level totals among courtship-access convergence cells ===")
    type_rows = sorted(by_type.items(), key=lambda kv: -(kv[1]["from_male"] * (kv[1]["from_PAM"] + kv[1]["from_OA"])))
    print(f"{'type':20s} {'n':>3s} {'Σm':>8s} {'Σf':>8s} {'ΣPAM':>8s} {'ΣOA':>8s} {'Σ→P1':>8s} {'Σ→court':>8s}")
    for t, d in type_rows[:30]:
        print(
            f"{t[:20]:20s} {type_n[t]:3d} {d['from_male']:8.0f} {d['from_female']:8.0f} "
            f"{d['from_PAM']:8.0f} {d['from_OA']:8.0f} {d['to_P1']:8.0f} {d['to_court']:8.0f}"
        )

    # Direct hop-1 male → courtship cells that ALSO receive mod
    print("\n=== courtship-core cells with direct male hop-1 input AND mod input ===")
    core_rows = []
    for t in court:
        mw = male_to.get(int(t), 0.0)
        if mw <= 0:
            continue
        pres, pw = ins(int(t))
        pam_w = oa_w = ppl_w = 0.0
        for p, ww in zip(pres, pw):
            p = int(p)
            tp = ntype(p)
            if tp.startswith("PAM"):
                pam_w += float(ww)
            elif tp.startswith("PPL"):
                ppl_w += float(ww)
            elif tp.startswith("OA-") or tp.startswith("OAN"):
                oa_w += float(ww)
        if pam_w + oa_w + ppl_w <= 0:
            continue
        core_rows.append(
            {
                "type": ntype(int(t)),
                "id": int(ids[t]),
                "from_male": mw,
                "from_female": female_to.get(int(t), 0.0),
                "from_PAM": pam_w,
                "from_OA": oa_w,
                "from_PPL": ppl_w,
            }
        )
    core_rows.sort(key=lambda r: -(r["from_male"] * (r["from_PAM"] + r["from_OA"] + r["from_PPL"])))
    print(f"n={len(core_rows)}")
    for r in core_rows[:30]:
        print(
            f"{r['type'][:18]:18s} id={r['id']} male={r['from_male']:.0f} "
            f"fem={r['from_female']:.0f} PAM={r['from_PAM']:.0f} OA={r['from_OA']:.0f} PPL={r['from_PPL']:.0f}"
        )

    # OA-VPM3 special: known MBON→P1 intermediate; check male + OA
    print("\n=== named intermediates of interest ===")
    for name in [
        "OA-VPM3",
        "OA-VUMa6",
        "SIP106m",
        "SMP163",
        "AOTU100m",
        "FLA001m",
        "SIP105m",
        "AN05B035",
        "IN05B011a",
        "AN09B017c",
        "AN09B017f",
        "AN03A008",
        "pC1_1a",
    ]:
        idxs = [i for i in range(len(ids)) if ntype(i) == name]
        if not idxs:
            print(name, "ABSENT")
            continue
        mw = sum(male_to.get(i, 0.0) for i in idxs)
        fw = sum(female_to.get(i, 0.0) for i in idxs)
        pam_w = oa_w = ppl_w = to_p1 = 0.0
        for i in idxs:
            pres, pw = ins(i)
            for p, ww in zip(pres, pw):
                p = int(p)
                tp = ntype(p)
                if tp.startswith("PAM"):
                    pam_w += float(ww)
                elif tp.startswith("PPL"):
                    ppl_w += float(ww)
                elif tp.startswith("OA-") or tp.startswith("OAN"):
                    oa_w += float(ww)
            tg, w = outs(i)
            for tt, ww in zip(tg, w):
                if int(tt) in p1set:
                    to_p1 += float(ww)
        print(
            f"{name:14s} n={len(idxs)} male={mw:.0f} fem={fw:.0f} "
            f"PAM={pam_w:.0f} OA={oa_w:.0f} PPL={ppl_w:.0f} →P1={to_p1:.0f}"
        )

    # Save compact JSON
    out = {
        "populations": {
            "male": int(len(male)),
            "female": int(len(female)),
            "PAM": int(len(pam)),
            "PPL": int(len(ppl)),
            "OA": int(len(oa)),
            "courtship_core": int(len(court)),
            "hop1_overlap_male_mod": int(len(both)),
            "overlap_with_court_access": int(len(court_acc)),
        },
        "top_convergence_cells": rows[:80],
        "type_level_court_access": [
            {"type": t, "n": type_n[t], **{k: float(v) for k, v in d.items()}}
            for t, d in type_rows[:40]
        ],
        "courtship_core_with_male_and_mod": core_rows[:80],
    }
    Path("convergence-report.json").write_text(json.dumps(out, indent=2, default=str))
    print("\nWrote convergence-report.json")


if __name__ == "__main__":
    main()
