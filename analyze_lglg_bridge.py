"""Focused bridge analysis: why LgLG5 > LgLG6 effective drive to P1."""
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


def main() -> None:
    meta = json.loads((ROOT / "full-meta.json").read_text())
    nodes = np.load(ROOT / "full-nodes.npz")
    ids, signs = nodes["ids"], nodes["signs"].copy()
    ann = feather.read_table(ROOT / "annotations.feather").to_pylist()
    by_id = {r["bodyId"]: r for r in ann}
    W = load_npz(ROOT / "full-graph.npz").tocsr()
    indptr, indices, data = W.indptr, W.indices, W.data
    p1 = np.searchsorted(ids, P1_IDS)
    p1set = set(int(x) for x in p1)
    vab3 = np.searchsorted(ids, VAB3_IDS)
    signs[vab3] = 1
    signs[np.asarray(meta["groups"]["candidate_female"], np.int32)] = 1

    def ntype(i):
        return (by_id.get(int(ids[i]), {}) or {}).get("type") or "?"

    def outs(i):
        s, e = indptr[i], indptr[i + 1]
        return indices[s:e], data[s:e]

    def cells(t):
        return np.asarray([i for i in range(len(ids)) if ntype(i) == t], dtype=np.int32)

    lg6, lg5 = cells("LgLG6"), cells("LgLG5")
    mal = np.asarray(meta["groups"]["mAL"], np.int32)

    # Types of interest: AN09B017*, shared hubs, P1 excitatory sources, mAL
    focus_types = [
        "AN09B017a", "AN09B017b", "AN09B017c", "AN09B017d",
        "AN09B017e", "AN09B017f", "AN09B017g",
        "AN05B035", "IN05B011a", "IN05B011b",
        "AN05B023a", "AN05B023b", "AN05B023c", "AN05B023d",
        "FLA001m", "SIP105m", "SIP025", "VES206m", "AVLP720m", "AVLP721m",
        "AN08B020", "pC1_1a",
        "mAL_m8", "mAL_m1", "mAL_m5b", "SIP100m", "SIP112m",
    ]

    def weight_between(src_pop, dst_pop):
        dset = set(int(x) for x in dst_pop)
        tot = 0.0
        n = 0
        for i in src_pop:
            tg, w = outs(int(i))
            for t, ww in zip(tg, w):
                if int(t) in dset:
                    tot += float(ww)
                    n += 1
        return tot, n

    def incoming_from_all_to_p1():
        """For each focus type: synapses to P1 and sign."""
        acc = {}
        for typ in focus_types:
            pop = cells(typ)
            if pop.size == 0:
                acc[typ] = None
                continue
            tot, n = weight_between(pop, p1)
            acc[typ] = {
                "n_cells": int(pop.size),
                "sign": int(signs[pop[0]]) if pop.size else None,
                "to_P1_syn": tot,
                "to_P1_edges": n,
                "from_LgLG6_syn": weight_between(lg6, pop)[0],
                "from_LgLG6_edges": weight_between(lg6, pop)[1],
                "from_LgLG5_syn": weight_between(lg5, pop)[0],
                "from_LgLG5_edges": weight_between(lg5, pop)[1],
            }
        return acc

    bridge = incoming_from_all_to_p1()
    print("=== focus types: source→type→P1 ===")
    print(f"{'type':16s} {'sg':>3s} {'n':>3s} {'6→T':>8s} {'5→T':>8s} {'T→P1':>8s} {'Δ6-5':>8s}")
    for typ, info in bridge.items():
        if not info:
            continue
        d = info["from_LgLG6_syn"] - info["from_LgLG5_syn"]
        print(
            f"{typ:16s} {info['sign']:+3d} {info['n_cells']:3d} "
            f"{info['from_LgLG6_syn']:8.0f} {info['from_LgLG5_syn']:8.0f} "
            f"{info['to_P1_syn']:8.0f} {d:+8.0f}"
        )

    # Per-cell average drive (normalize by source population size)
    print("\n=== per-source-cell mean synapses into focus types ===")
    print(f"{'type':16s} {'6/cell':>8s} {'5/cell':>8s} {'ratio 6/5':>10s}")
    for typ, info in bridge.items():
        if not info:
            continue
        a = info["from_LgLG6_syn"] / len(lg6)
        b = info["from_LgLG5_syn"] / len(lg5)
        ratio = a / b if b > 0 else float("inf")
        print(f"{typ:16s} {a:8.1f} {b:8.1f} {ratio:10.3f}")

    # mAL access
    print("\n=== access to mAL (inhibitory P1 brake) ===")
    print("LgLG6→mAL", weight_between(lg6, mal))
    print("LgLG5→mAL", weight_between(lg5, mal))
    print("per-cell 6", weight_between(lg6, mal)[0] / len(lg6), "5", weight_between(lg5, mal)[0] / len(lg5))

    # Who among hop-1 of each provides most EXCITATORY product to P1
    # (source +, intermediate sign s, path to P1)
    print("\n=== hop-1 intermediates ranked by signed 2-hop product to P1 ===")

    def rank_intermediates(src_pop, label):
        rows = []
        for i in src_pop:
            tg, w = outs(int(i))
            for t, ww in zip(tg, w):
                t = int(t)
                if signs[int(i)] == 0 or signs[t] == 0:
                    continue
                # direct t→P1
                tg2, w2 = outs(t)
                for t2, ww2 in zip(tg2, w2):
                    if int(t2) in p1set:
                        prod = float(ww) * float(ww2) * np.sign(signs[int(i)]) * np.sign(signs[t])
                        rows.append((ntype(t), t, prod, float(ww), float(ww2), "d2"))
        by = defaultdict(lambda: {"prod": 0.0, "n": 0, "in_w": 0.0, "out_w": 0.0})
        for typ, _, prod, w1, w2, _ in rows:
            by[typ]["prod"] += prod
            by[typ]["n"] += 1
            by[typ]["in_w"] += w1
            by[typ]["out_w"] += w2
        print(f"\n{label} 2-hop via intermediate → P1 (top 12 by |prod|):")
        for typ, d in sorted(by.items(), key=lambda kv: -abs(kv[1]["prod"]))[:12]:
            print(
                f"  {typ:16s} prod={d['prod']:12.0f} n={d['n']:4d} "
                f"in_syn={d['in_w']:8.0f} out_syn={d['out_w']:8.0f} "
                f"sign={int(signs[cells(typ)[0]]) if cells(typ).size else '?'}"
            )
        return by

    r6 = rank_intermediates(lg6, "LgLG6")
    r5 = rank_intermediates(lg5, "LgLG5")

    print("\n=== 2-hop product ratio LgLG6/LgLG5 by type (only both nonzero) ===")
    all_types = sorted(set(r6) | set(r5))
    for typ in all_types:
        p6 = r6.get(typ, {}).get("prod", 0.0)
        p5 = r5.get(typ, {}).get("prod", 0.0)
        if p5 == 0 and p6 == 0:
            continue
        ratio = p6 / p5 if p5 != 0 else float("inf")
        print(f"  {typ:16s} p6={p6:12.0f} p5={p5:12.0f} ratio={ratio:8.3f}")

    # Cross-inhibition between LgLG6 and LgLG5
    print("\n=== LgLG6 ↔ LgLG5 crosstalk ===")
    print("6→5", weight_between(lg6, lg5), "5→6", weight_between(lg5, lg6))
    print(
        "per-cell 6→5",
        weight_between(lg6, lg5)[0] / len(lg6),
        "5→6",
        weight_between(lg5, lg6)[0] / len(lg5),
    )

    # Sum of excitatory vs inhibitory P1 source weight reachable at hop2
    print("\n=== hop-2 mass into P1-excitatory vs P1-inhibitory sources ===")
    # classify all neurons by sign and whether they touch P1
    p1_sources_e = []
    p1_sources_i = []
    Wcsc = W.tocsc()
    for t in p1:
        s, e = Wcsc.indptr[t], Wcsc.indptr[t + 1]
        for k in range(s, e):
            p = int(Wcsc.indices[k])
            w = float(Wcsc.data[k])
            if signs[p] > 0:
                p1_sources_e.append((p, w))
            elif signs[p] < 0:
                p1_sources_i.append((p, w))
    e_set = set(p for p, _ in p1_sources_e)
    i_set = set(p for p, _ in p1_sources_i)
    e_wmap = defaultdict(float)
    i_wmap = defaultdict(float)
    for p, w in p1_sources_e:
        e_wmap[p] += w
    for p, w in p1_sources_i:
        i_wmap[p] += w

    def hop2_into(src_pop, target_wmap):
        # sum over src→t of w(src,t) * target_wmap[t] for t in target_wmap
        tot = 0.0
        by_type = defaultdict(float)
        for i in src_pop:
            tg, w = outs(int(i))
            for t, ww in zip(tg, w):
                t = int(t)
                if t in target_wmap:
                    contrib = float(ww) * target_wmap[t]
                    tot += contrib
                    by_type[ntype(t)] += contrib
        return tot, sorted(by_type.items(), key=lambda kv: -kv[1])[:10]

    for label, pop in (("LgLG6", lg6), ("LgLG5", lg5)):
        te, top_e = hop2_into(pop, e_wmap)
        ti, top_i = hop2_into(pop, i_wmap)
        print(f"\n{label} hop2 into P1-exc sources: {te:.0f} top={top_e}")
        print(f"{label} hop2 into P1-inh sources: {ti:.0f} top={top_i}")
        print(f"{label} net (exc-inh) P1-source weighted access: {te - ti:.0f}")

    out = {
        "bridge": bridge,
        "LgLG6_2hop": {k: dict(v) for k, v in r6.items()},
        "LgLG5_2hop": {k: dict(v) for k, v in r5.items()},
    }
    Path("pathway-bridge-lg6-lg5.json").write_text(json.dumps(out, indent=2, default=str))
    print("\nWrote pathway-bridge-lg6-lg5.json")


if __name__ == "__main__":
    main()
