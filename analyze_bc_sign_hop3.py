"""AN09B017b/c sign evidence + LgLG5/6 hop-3 path decomposition.

Static anatomical analysis only. No outcome-seeking simulation.
Outputs: bc-hop3-evidence.json
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
VAB3_IDS = [11998, 13341, 13693, 512498]
BC_TYPES = ("AN09B017b", "AN09B017c")


def main() -> None:
    meta = json.loads((ROOT / "full-meta.json").read_text())
    nodes = np.load(ROOT / "full-nodes.npz")
    ids, base_signs = nodes["ids"], nodes["signs"].copy()
    ann = feather.read_table(ROOT / "annotations.feather").to_pylist()
    by_id = {r["bodyId"]: r for r in ann}
    nt = feather.read_table(ROOT / "nt.feather").to_pydict()
    ntmap = {b: i for i, b in enumerate(nt["body"])}
    W = load_npz(ROOT / "full-graph.npz").tocsr()
    indptr, indices, data = W.indptr, W.indices, W.data

    def ntype(i: int) -> str:
        return (by_id.get(int(ids[i]), {}) or {}).get("type") or "?"

    def outs(i: int):
        s, e = indptr[i], indptr[i + 1]
        return indices[s:e], data[s:e]

    def cells(t: str) -> np.ndarray:
        return np.asarray([i for i in range(len(ids)) if ntype(i) == t], dtype=np.int32)

    p1 = np.searchsorted(ids, P1_IDS)
    p1set = set(int(x) for x in p1)
    vab3 = np.searchsorted(ids, VAB3_IDS)

    def apply_signs(bc_sign: int) -> np.ndarray:
        """Baseline overrides + optional AN09B017b/c sign."""
        signs = base_signs.copy()
        signs[vab3] = 1
        signs[np.asarray(meta["groups"]["candidate_female"], np.int32)] = 1
        for t in BC_TYPES:
            idx = cells(t)
            signs[idx] = bc_sign
        return signs

    # --- evidence ledger for AN09B017* ---
    ledger = {}
    for t in [
        "AN09B017a",
        "AN09B017b",
        "AN09B017c",
        "AN09B017d",
        "AN09B017e",
        "AN09B017f",
        "AN09B017g",
    ]:
        rows = [r for r in ann if r.get("type") == t]
        entries = []
        for r in rows:
            bid = r["bodyId"]
            i = ntmap.get(bid)
            entries.append(
                {
                    "bodyId": bid,
                    "dimorphism": r.get("dimorphism"),
                    "fruDsx": r.get("fruDsx"),
                    "superclass": r.get("superclass"),
                    "synonyms": r.get("synonyms") or "",
                    "receptorType": r.get("receptorType"),
                    "hemibrainType": r.get("hemibrainType"),
                    "mancType": r.get("mancType"),
                    "consensus_nt": nt["consensus_nt"][i] if i is not None else None,
                    "predicted_nt": nt["predicted_nt"][i] if i is not None else None,
                    "predicted_nt_confidence": (
                        float(nt["predicted_nt_confidence"][i]) if i is not None else None
                    ),
                    "celltype_predicted_nt": (
                        nt["celltype_predicted_nt"][i] if i is not None else None
                    ),
                    "celltype_predicted_nt_confidence": (
                        float(nt["celltype_predicted_nt_confidence"][i])
                        if i is not None
                        else None
                    ),
                    "ground_truth_nt": nt["ground_truth"][i] if i is not None else None,
                }
            )
        idx = cells(t)
        model_sign = int(base_signs[idx[0]]) if idx.size else None
        after_override = int(apply_signs(-1)[idx[0]]) if idx.size else None
        to_p1 = 0.0
        for i in idx:
            tg, w = outs(int(i))
            for tt, ww in zip(tg, w):
                if int(tt) in p1set:
                    to_p1 += float(ww)
        ledger[t] = {
            "n_cells": int(idx.size),
            "bodyIds": [int(ids[i]) for i in idx],
            "model_default_sign": model_sign,
            "baseline_sign_after_vAB3_override": after_override,
            "to_P1_syn": to_p1,
            "cells": entries,
        }

    # claim classification
    claim = {
        "anatomical_fact": [
            "AN09B017b/c are male-specific fru_high ascending neurons (2 cells each).",
            "consensus_nt and ground_truth are glutamate for all AN09B017 subtypes.",
            "receptorType is empty for every AN09B017 cell — no postsynaptic receptor annotation in this dataset.",
            "LgLG6→AN09B017b/c synapses ≫ LgLG5 (prior hop-2 analysis).",
            "Direct b/c→P1 synapses are small (b: 7, c: 55) relative to fan-out to mAL / AVLP / SIP.",
        ],
        "model_assumption": [
            "Fast sign map: glutamate → −1 (inhibitory). This is why b/c act as a brake.",
            "vAB3e/f (synonym Yu 2010) are forced +1 despite the same consensus glutamate.",
            "AN09B017g carries literature vAB3 synonyms but is NOT forced +1 in this model.",
            "Predicted_nt for b/c is acetylcholine at conf≈0.94–0.96, disagreeing with consensus glutamate; the model uses consensus, not predicted.",
            "Uniform Poisson sensory encoding on LgLG6/7 vs LgLG5/8.",
        ],
        "biological_hypothesis": [
            "Glutamate release from b/c is inhibitory onto P1-related targets in vivo (would require GluCl-type postsynaptic receptors).",
            "Predicted ACh would make b/c excitatory and reverse their circuit role.",
            "b/c loading asymmetry between LgLG6 and LgLG5 is a functional male-cue brake, not just a wiring footprint.",
        ],
        "not_established": [
            "That AN09B017b/c inhibit P1 in the animal.",
            "That relieving b/c changes mate preference.",
            "That consensus glutamate is the functional sign at these synapses.",
        ],
    }

    # --- hop-2 / hop-3 signed products under sign variants ---
    lg6, lg5 = cells("LgLG6"), cells("LgLG5")

    def hop2_p1_source_access(src_pop, signs):
        """Weighted access into P1-excitatory vs P1-inhibitory sources."""
        p1_sources_e = defaultdict(float)
        p1_sources_i = defaultdict(float)
        Wcsc = W.tocsc()
        for t in p1:
            s, e = Wcsc.indptr[t], Wcsc.indptr[t + 1]
            for k in range(s, e):
                p = int(Wcsc.indices[k])
                w = float(Wcsc.data[k])
                if signs[p] > 0:
                    p1_sources_e[p] += w
                elif signs[p] < 0:
                    p1_sources_i[p] += w
        te = ti = 0.0
        by_i = defaultdict(float)
        for i in src_pop:
            tg, w = outs(int(i))
            for t, ww in zip(tg, w):
                t = int(t)
                if t in p1_sources_e:
                    te += float(ww) * p1_sources_e[t]
                if t in p1_sources_i:
                    c = float(ww) * p1_sources_i[t]
                    ti += c
                    by_i[ntype(t)] += c
        return {
            "into_P1_exc_sources": te,
            "into_P1_inh_sources": ti,
            "net": te - ti,
            "top_inh_types": sorted(by_i.items(), key=lambda kv: -kv[1])[:8],
        }

    def hop3_decomposition(src_pop, signs, label):
        """Attribute depth-3 signed product mass to (t1_type, t2_type).

        Paths: src → t1 → t2 → P1
        product = w1 * w2 * w3 * sign(src)*sign(t1)*sign(t2)
        """
        by_t1 = defaultdict(float)
        by_t2 = defaultdict(float)
        by_pair = defaultdict(float)
        pos = neg = 0.0
        n_paths = 0
        for i in src_pop:
            sg0 = int(signs[int(i)])
            if sg0 == 0:
                continue
            tg1, w1 = outs(int(i))
            for t1, ww1 in zip(tg1, w1):
                t1 = int(t1)
                sg1 = int(signs[t1])
                if sg1 == 0:
                    continue
                tg2, w2 = outs(t1)
                for t2, ww2 in zip(tg2, w2):
                    t2 = int(t2)
                    sg2 = int(signs[t2])
                    if sg2 == 0:
                        continue
                    tg3, w3 = outs(t2)
                    for t3, ww3 in zip(tg3, w3):
                        if int(t3) not in p1set:
                            continue
                        prod = (
                            float(ww1)
                            * float(ww2)
                            * float(ww3)
                            * sg0
                            * sg1
                            * sg2
                        )
                        n_paths += 1
                        if prod > 0:
                            pos += prod
                        else:
                            neg += prod
                        ty1, ty2 = ntype(t1), ntype(t2)
                        by_t1[ty1] += prod
                        by_t2[ty2] += prod
                        by_pair[(ty1, ty2)] += prod
        return {
            "label": label,
            "n_paths": n_paths,
            "pos_mass": pos,
            "neg_mass": neg,
            "net_mass": pos + neg,
            "by_t1": sorted(by_t1.items(), key=lambda kv: -kv[1]),
            "by_t2": sorted(by_t2.items(), key=lambda kv: -kv[1]),
            "by_pair_top": sorted(by_pair.items(), key=lambda kv: -abs(kv[1]))[:30],
        }

    sign_variants = {
        "bc_inh_minus1": apply_signs(-1),
        "bc_exc_plus1": apply_signs(1),
        "bc_zero": apply_signs(0),
    }

    access = {}
    hop3 = {}
    for vname, signs in sign_variants.items():
        access[vname] = {
            "LgLG6": hop2_p1_source_access(lg6, signs),
            "LgLG5": hop2_p1_source_access(lg5, signs),
        }
        hop3[vname] = {
            "LgLG6": hop3_decomposition(lg6, signs, "LgLG6"),
            "LgLG5": hop3_decomposition(lg5, signs, "LgLG5"),
        }

    # rank female-specific amplifiers / male sinks / shared under baseline signs
    base = sign_variants["bc_inh_minus1"]
    d6 = hop3["bc_inh_minus1"]["LgLG6"]
    d5 = hop3["bc_inh_minus1"]["LgLG5"]
    m6 = dict(d6["by_t1"])
    m5 = dict(d5["by_t1"])
    all_t1 = set(m6) | set(m5)
    ranks = []
    for t in all_t1:
        p6, p5 = m6.get(t, 0.0), m5.get(t, 0.0)
        ranks.append(
            {
                "t1_type": t,
                "LgLG6": p6,
                "LgLG5": p5,
                "delta_5_minus_6": p5 - p6,
                "abs_ratio_5_over_6": (abs(p5) / abs(p6)) if p6 else None,
            }
        )
    ranks_sorted = sorted(ranks, key=lambda r: -r["delta_5_minus_6"])
    male_sink = sorted(ranks, key=lambda r: r["delta_5_minus_6"])  # most male-negative or female-weaker
    shared = [
        r
        for r in ranks
        if abs(r["LgLG6"]) > 1e3 and abs(r["LgLG5"]) > 1e3 and abs(r["delta_5_minus_6"]) < 0.5 * (abs(r["LgLG6"]) + abs(r["LgLG5"]))
    ]
    shared = sorted(shared, key=lambda r: -(abs(r["LgLG6"]) + abs(r["LgLG5"])))[:15]

    # central switch candidates: types that (a) receive hop-1 mass asymmetrically
    # and (b) have large |signed mass| as t1 into P1 paths, and (c) touch P1-relevant structure
    # Also compute hop-1 syn into each focus type for LgLG6 vs LgLG5
    focus = [
        "AN09B017a",
        "AN09B017b",
        "AN09B017c",
        "AN09B017d",
        "AN09B017e",
        "AN09B017f",
        "AN09B017g",
        "AN05B035",
        "IN05B011a",
        "IN05B011b",
        "FLA001m",
        "SIP105m",
        "SIP025",
        "VES206m",
        "AVLP720m",
        "AVLP721m",
        "AN08B020",
        "AN03A008",
        "pC1_1a",
        "mAL_m1",
        "mAL_m8",
        "SIP100m",
        "SIP112m",
        "SIP113m",
        "VES022",
        "LH004m",
        "PVLP048",
        "GNG700m",
    ]

    def syn_between(src, dst):
        dset = set(int(x) for x in dst)
        tot = n = 0.0
        for i in src:
            tg, w = outs(int(i))
            for t, ww in zip(tg, w):
                if int(t) in dset:
                    tot += float(ww)
                    n += 1
        return tot, n

    switch_rows = []
    for typ in focus:
        pop = cells(typ)
        if pop.size == 0:
            continue
        w6, e6 = syn_between(lg6, pop)
        w5, e5 = syn_between(lg5, pop)
        # direct to P1
        to_p1 = 0.0
        for i in pop:
            tg, w = outs(int(i))
            for t, ww in zip(tg, w):
                if int(t) in p1set:
                    to_p1 += float(ww)
        # hop-1 product mass as t1 (from hop3 by_t1)
        p6 = m6.get(typ, 0.0)
        p5 = m5.get(typ, 0.0)
        switch_rows.append(
            {
                "type": typ,
                "sign": int(base[pop[0]]),
                "n": int(pop.size),
                "from_LgLG6_syn": w6,
                "from_LgLG5_syn": w5,
                "ratio_6_over_5": (w6 / w5) if w5 else None,
                "to_P1_syn": to_p1,
                "hop3_mass_as_t1_LgLG6": p6,
                "hop3_mass_as_t1_LgLG5": p5,
            }
        )

    report = {
        "evidence_ledger": ledger,
        "claim_classes": claim,
        "hop2_P1_source_access_by_sign_variant": access,
        "hop3_by_sign_variant": {
            k: {
                "LgLG6": {kk: vv for kk, vv in v["LgLG6"].items() if kk != "by_pair_top"},
                "LgLG5": {kk: vv for kk, vv in v["LgLG5"].items() if kk != "by_pair_top"},
                "LgLG6_pair_top": v["LgLG6"]["by_pair_top"][:20],
                "LgLG5_pair_top": v["LgLG5"]["by_pair_top"][:20],
            }
            for k, v in hop3.items()
        },
        "female_amplifier_rank": ranks_sorted[:20],
        "male_sink_rank": male_sink[:20],
        "shared_relay_rank": shared,
        "central_switch_candidates": sorted(
            switch_rows, key=lambda r: -abs(r["hop3_mass_as_t1_LgLG5"] - r["hop3_mass_as_t1_LgLG6"])
        ),
    }
    Path("bc-hop3-evidence.json").write_text(json.dumps(report, indent=2, default=str))

    # console summary
    print("=== AN09B017 sign evidence (b/c focus) ===")
    for t in BC_TYPES + ("AN09B017e", "AN09B017f", "AN09B017g", "AN09B017d"):
        L = ledger[t]
        c0 = L["cells"][0]
        print(
            f"{t}: model={L['model_default_sign']} consensus={c0['consensus_nt']} "
            f"pred={c0['predicted_nt']}({c0['predicted_nt_confidence']}) "
            f"gt={c0['ground_truth_nt']} receptor={c0['receptorType']} "
            f"syn={c0['synonyms']!r} →P1={L['to_P1_syn']}"
        )

    print("\n=== hop-2 P1-source access by bc sign ===")
    for v, d in access.items():
        print(v, "LgLG6 net", round(d["LgLG6"]["net"]), "LgLG5 net", round(d["LgLG5"]["net"]))

    print("\n=== hop-3 signed product mass by bc sign ===")
    for v, d in hop3.items():
        print(
            v,
            "LgLG6",
            {k: d["LgLG6"][k] for k in ("n_paths", "pos_mass", "neg_mass", "net_mass")},
            "LgLG5",
            {k: d["LgLG5"][k] for k in ("n_paths", "pos_mass", "neg_mass", "net_mass")},
        )

    print("\n=== top female amplifiers (t1, delta 5-6) ===")
    for r in ranks_sorted[:12]:
        print(f"  {r['t1_type']:16s} 6={r['LgLG6']:12.0f} 5={r['LgLG5']:12.0f} d={r['delta_5_minus_6']:12.0f}")

    print("\n=== top male sinks / male-heavier t1 ===")
    for r in male_sink[:8]:
        print(f"  {r['t1_type']:16s} 6={r['LgLG6']:12.0f} 5={r['LgLG5']:12.0f} d={r['delta_5_minus_6']:12.0f}")

    print("\n=== central switch candidates (|Δ hop3 mass|) ===")
    for r in report["central_switch_candidates"][:12]:
        print(
            f"  {r['type']:16s} sg={r['sign']:+d} 6→T={r['from_LgLG6_syn']:8.0f} "
            f"5→T={r['from_LgLG5_syn']:8.0f} →P1={r['to_P1_syn']:6.0f} "
            f"hop3_6={r['hop3_mass_as_t1_LgLG6']:12.0f} hop3_5={r['hop3_mass_as_t1_LgLG5']:12.0f}"
        )

    print("\nWrote bc-hop3-evidence.json")


if __name__ == "__main__":
    main()
