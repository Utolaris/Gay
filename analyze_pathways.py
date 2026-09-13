"""Static signed pathway audit on the prepared MaleCNS graph (read-only)."""
from __future__ import annotations

import json
from collections import Counter, defaultdict
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
    ids, base_signs = nodes["ids"], nodes["signs"]
    W = load_npz(ROOT / "full-graph.npz").tocsr()
    ann = feather.read_table(ROOT / "annotations.feather").to_pylist()
    by_id = {r["bodyId"]: r for r in ann}

    signs = base_signs.copy()
    p1 = np.searchsorted(ids, P1_IDS)
    vab3 = np.searchsorted(ids, VAB3_IDS)
    mal = np.asarray(meta["groups"]["mAL"], np.int32)
    male = np.asarray(meta["groups"]["candidate_male"], np.int32)
    female = np.asarray(meta["groups"]["candidate_female"], np.int32)
    broad = np.asarray(meta["groups"]["P1_related"], np.int32)
    signs[vab3] = 1
    signs[female] = 1  # baseline override under audit

    def ntype(i: int) -> str:
        return (by_id.get(int(ids[i]), {}) or {}).get("type") or "?"

    def nts(i: int) -> str:
        return (by_id.get(int(ids[i]), {}) or {}).get("consensus_nt") or "?"

    # annotations may not carry nt; use signs + type
    indptr, indices, data = W.indptr, W.indices, W.data

    def outs(i: int):
        s, e = indptr[i], indptr[i + 1]
        return indices[s:e], data[s:e]

    def incoming_to(targets: np.ndarray) -> dict:
        acc = defaultdict(float)
        edges = defaultdict(int)
        tset = set(int(x) for x in targets)
        # iterate all neurons' outgoing — expensive but graph is 25M; use CSC-like via transpose once
        return acc, edges

    # Build reverse adjacency only for P1 column via csc
    Wcsc = W.tocsc()
    print("=== direct presynaptic classes -> P1 (weighted, nonzero fast sign) ===")
    acc = defaultdict(float)
    edge_n = defaultdict(int)
    for t in p1:
        start, end = Wcsc.indptr[t], Wcsc.indptr[t + 1]
        for k in range(start, end):
            p = Wcsc.indices[k]
            w = Wcsc.data[k]
            if signs[p] == 0:
                continue
            key = (ntype(p), int(signs[p]))
            acc[key] += float(w)
            edge_n[key] += 1
    for (typ, sg), w in sorted(acc.items(), key=lambda kv: -kv[1])[:30]:
        print(f"  {typ:28s} sign={sg:+d} edges={edge_n[(typ, sg)]:4d} syn={w:.0f}")

    def out_by_type(pre: np.ndarray, label: str) -> None:
        print(f"\n=== 1-hop out from {label} (n={len(pre)}) ===")
        a = defaultdict(float)
        en = defaultdict(int)
        for i in pre:
            if signs[i] == 0:
                continue
            tg, w = outs(i)
            for t, ww in zip(tg, w):
                key = (ntype(t), int(signs[i]))
                a[key] += float(ww)
                en[key] += 1
        for (typ, sg), w in sorted(a.items(), key=lambda kv: -kv[1])[:20]:
            print(f"  -> {typ:28s} pre_sign={sg:+d} edges={en[(typ, sg)]:4d} syn={w:.0f}")

    out_by_type(male, "candidate_male")
    out_by_type(female, "candidate_female")

    p1set = set(int(x) for x in p1)
    malset = set(int(x) for x in mal)

    def path_stats(pre: np.ndarray, label: str) -> None:
        direct = 0.0
        d_edges = 0
        to_mal = 0.0
        tm_edges = 0
        depth2 = defaultdict(float)  # (j_type, j_sign, path_sign) -> product weight
        for i in pre:
            tg, w = outs(i)
            for j, w1 in zip(tg, w):
                if signs[i] == 0:
                    continue
                if j in p1set:
                    direct += float(w1)
                    d_edges += 1
                if j in malset:
                    to_mal += float(w1)
                    tm_edges += 1
                if signs[j] == 0:
                    continue
                tg2, w2 = outs(j)
                for t, wj in zip(tg2, w2):
                    if t in p1set:
                        key = (ntype(j), int(signs[j]), int(signs[i]) * int(signs[j]))
                        depth2[key] += float(w1) * float(wj)
        print(f"\n=== {label} -> structure ===")
        print(f"  direct -> P1: edges={d_edges} synapses={direct}")
        print(f"  direct -> mAL: edges={tm_edges} synapses={to_mal}")
        print("  depth-2 via intermediate (product of anatomical counts; not physiology):")
        for k, v in sorted(depth2.items(), key=lambda kv: -kv[1])[:15]:
            print(f"    via {k[0]:28s} j_sign={k[1]:+d} path_sign={k[2]:+d}  prod={v:.1f}")

    path_stats(male, "candidate_male")
    path_stats(female, "candidate_female")
    path_stats(vab3, "vAB3")

    print("\n=== mAL -> P1 direct ===")
    tot = 0.0
    n = 0
    for i in mal:
        tg, w = outs(i)
        for t, ww in zip(tg, w):
            if t in p1set:
                tot += float(ww)
                n += 1
    print(f"  edges={n} synapses={tot} mAL_sign_set={sorted(set(int(x) for x in signs[mal]))}")

    print("\n=== input population type composition ===")
    for label, arr in [
        ("male", male),
        ("female", female),
        ("mAL", mal),
        ("P1", p1),
        ("vAB3", vab3),
        ("broad_P1", broad),
    ]:
        c = Counter(ntype(int(i)) for i in arr)
        print(f"  {label}: {dict(c)}")

    # who provides most inhibition to P1
    print("\n=== inhibitory (sign=-1) sources into P1, top types ===")
    inh = defaultdict(float)
    for t in p1:
        start, end = Wcsc.indptr[t], Wcsc.indptr[t + 1]
        for k in range(start, end):
            p = Wcsc.indices[k]
            w = Wcsc.data[k]
            if signs[p] == -1:
                inh[ntype(p)] += float(w)
    for typ, w in sorted(inh.items(), key=lambda kv: -kv[1])[:15]:
        print(f"  {typ:28s} {w:.0f}")

    # bottleneck: high-outdegree neurons that both receive from male and touch P1 within 2 hops
    print("\n=== candidate bottleneck types (male 1-hop ∩ P1 1-in or 2-hop) ===")
    male_targets = set()
    for i in male:
        tg, _ = outs(i)
        male_targets.update(int(x) for x in tg)
    # P1 sources 1-hop
    p1_sources = set()
    for t in p1:
        start, end = Wcsc.indptr[t], Wcsc.indptr[t + 1]
        p1_sources.update(int(Wcsc.indices[k]) for k in range(start, end))
    bottleneck = male_targets & p1_sources
    bc = Counter(ntype(i) for i in bottleneck)
    print(f"  n_shared={len(bottleneck)} top types: {bc.most_common(20)}")
    # with sign product male->b->P1
    prod = defaultdict(float)
    for i in male:
        tg, w = outs(i)
        for j, w1 in zip(tg, w):
            if j not in bottleneck or signs[i] == 0 or signs[j] == 0:
                continue
            # weight j->P1
            tj, wj = outs(j)
            for t, w2 in zip(tj, wj):
                if t in p1set:
                    prod[(ntype(j), int(signs[i]) * int(signs[j]))] += float(w1) * float(w2)
    print("  product weights by (type, path_sign):")
    for k, v in sorted(prod.items(), key=lambda kv: -kv[1])[:15]:
        print(f"    {k[0]:28s} path_sign={k[1]:+d}  {v:.1f}")

    # Without female override, female fast weight
    print("\n=== female override impact ===")
    print("  base signs of female cells:", Counter(int(x) for x in base_signs[female]))
    print("  override signs:", Counter(int(x) for x in signs[female]))
    print("  male base signs:", Counter(int(x) for x in base_signs[male]))


if __name__ == "__main__":
    main()
