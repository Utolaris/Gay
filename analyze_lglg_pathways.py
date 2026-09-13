"""Signed pathway decomposition: LgLG6 (male) vs LgLG5 (female) → P1.

Anatomical only. No simulation. Reports per-hop excitatory/inhibitory mass,
edge counts, fan-out, recurrence, and the first hop where relative drive
diverges strongly.
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
MAX_HOP = 4
# Cap frontier expansion for hop>=3 (high-degree hubs explode).
FRONTIER_CAP = 5000


def main() -> None:
    meta = json.loads((ROOT / "full-meta.json").read_text())
    nodes = np.load(ROOT / "full-nodes.npz")
    ids, base_signs = nodes["ids"], nodes["signs"]
    signs = base_signs.copy()
    ann = feather.read_table(ROOT / "annotations.feather").to_pylist()
    by_id = {r["bodyId"]: r for r in ann}
    W = load_npz(ROOT / "full-graph.npz").tocsr()
    indptr, indices, data = W.indptr, W.indices, W.data

    p1 = np.searchsorted(ids, P1_IDS)
    vab3 = np.searchsorted(ids, VAB3_IDS)
    signs[vab3] = 1
    # baseline female override as in sensory experiments
    female_all = np.asarray(meta["groups"]["candidate_female"], np.int32)
    signs[female_all] = 1

    def ntype(i: int) -> str:
        return (by_id.get(int(ids[i]), {}) or {}).get("type") or "?"

    def type_of_set(idxs) -> dict:
        c = defaultdict(int)
        for i in idxs:
            c[ntype(int(i))] += 1
        return dict(c)

    lg6 = np.asarray(
        [i for i in range(len(ids)) if ntype(i) == "LgLG6"], dtype=np.int32
    )
    lg5 = np.asarray(
        [i for i in range(len(ids)) if ntype(i) == "LgLG5"], dtype=np.int32
    )
    lg7 = np.asarray(
        [i for i in range(len(ids)) if ntype(i) == "LgLG7"], dtype=np.int32
    )
    lg8 = np.asarray(
        [i for i in range(len(ids)) if ntype(i) == "LgLG8"], dtype=np.int32
    )
    p1set = set(int(x) for x in p1)

    print("=== population sizes ===")
    for name, arr in [("LgLG6", lg6), ("LgLG5", lg5), ("LgLG7", lg7), ("LgLG8", lg8)]:
        print(name, len(arr), "signs", dict(zip(*np.unique(signs[arr], return_counts=True))))

    def outs(i: int):
        s, e = indptr[i], indptr[i + 1]
        return indices[s:e], data[s:e]

    def hop_stats(pre: np.ndarray, hop: int, cap: int | None = None) -> dict:
        """One-hop expansion: total E/I weight, edges, unique targets, fan-out."""
        e_w = i_w = 0.0
        e_n = i_n = 0
        targets = set()
        fanouts = []
        # if pre is large, sample or use all — LgLG6 is 16 cells, fine
        for i in pre:
            tg, w = outs(int(i))
            if tg.size == 0:
                fanouts.append(0)
                continue
            fanouts.append(int(tg.size))
            sg = signs[int(i)]
            if sg == 0:
                continue
            for t, ww in zip(tg, w):
                targets.add(int(t))
                if sg > 0:
                    e_w += float(ww)
                    e_n += 1
                else:
                    i_w += float(ww)
                    i_n += 1
        tlist = np.fromiter(targets, dtype=np.int32) if targets else np.zeros(0, np.int32)
        if cap is not None and tlist.size > cap:
            # keep highest in-degree-ish: prefer types that also touch P1 later
            # deterministic: sort by index, keep first cap — document truncation
            tlist = np.sort(tlist)[:cap]
        return {
            "hop": hop,
            "n_pre": int(len(pre)),
            "exc_syn": e_w,
            "inh_syn": i_w,
            "exc_edges": e_n,
            "inh_edges": i_n,
            "n_targets": int(len(targets)) if cap is None else int(tlist.size),
            "targets_truncated": cap is not None and len(targets) > cap,
            "mean_fanout": float(np.mean(fanouts)) if fanouts else 0.0,
            "median_fanout": float(np.median(fanouts)) if fanouts else 0.0,
            "max_fanout": int(np.max(fanouts)) if fanouts else 0,
            "target_set": tlist,
            "target_types_top": sorted(
                type_of_set(tlist).items(), key=lambda kv: -kv[1]
            )[:12],
        }

    def path_product_to_p1(pre: np.ndarray, max_depth: int) -> dict:
        """DFS-like layered product mass from pre cells to P1 (anatomical products)."""
        # layer 0: pre
        # accumulate product weights of paths that reach P1 at each depth
        hits_by_depth = defaultdict(float)  # depth -> sum of products
        path_counts = defaultdict(int)
        # BFS of (node, product_so_far, depth) with pruning
        # state: current frontier as dict node -> (sum of path products reaching it, depth)
        frontier = {int(i): 1.0 for i in pre}
        visited_products = defaultdict(float)  # node -> product mass (merged)
        for i in pre:
            visited_products[int(i)] += 1.0
        for depth in range(1, max_depth + 1):
            nxt = defaultdict(float)
            for node, prod in frontier.items():
                tg, w = outs(node)
                sg = signs[node]
                if sg == 0:
                    continue
                for t, ww in zip(tg, w):
                    p2 = prod * float(ww) * (1.0 if sg > 0 else -1.0)
                    t = int(t)
                    if t in p1set:
                        hits_by_depth[depth] += p2
                        path_counts[depth] += 1
                    nxt[t] += p2
            # merge into frontier (allow revisit for recurrence mass)
            frontier = dict(nxt)
            if len(frontier) > 200000:
                # keep largest mass nodes
                top = sorted(frontier.items(), key=lambda kv: -abs(kv[1]))[:50000]
                frontier = dict(top)
        return {
            "signed_product_by_depth": {int(k): v for k, v in hits_by_depth.items()},
            "path_edge_hits_by_depth": {int(k): v for k, v in path_counts.items()},
            "final_frontier_size": len(frontier),
        }

    def recurrence_stats(pre: np.ndarray, k: int = 3) -> dict:
        """Fraction of k-hop neighborhood that returns to the source set."""
        src = set(int(x) for x in pre)
        front = set(src)
        for _ in range(k):
            nxt = set()
            for i in front:
                tg, _ = outs(i)
                nxt.update(int(x) for x in tg)
            front = nxt
        return {
            "k": k,
            "neighborhood": len(front),
            "overlap_with_source": len(front & src),
            "source_in_k_hop_out": len(src & front) / max(len(src), 1),
        }

    def compare(label_a, a, label_b, b):
        print(f"\n######## {label_a} vs {label_b} ########")
        rows = {}
        for label, pop in ((label_a, a), (label_b, b)):
            hops = []
            cur = pop
            for h in range(1, MAX_HOP + 1):
                st = hop_stats(cur, h, cap=FRONTIER_CAP if h >= 3 else None)
                hops.append(st)
                cur = st["target_set"]
                if cur.size == 0:
                    break
            rows[label] = hops
            print(f"\n--- {label} n={len(pop)} ---")
            for st in hops:
                print(
                    f"  hop{st['hop']}: pre={st['n_pre']:5d} "
                    f"exc_syn={st['exc_syn']:10.0f} inh_syn={st['inh_syn']:10.0f} "
                    f"net={st['exc_syn']-st['inh_syn']:10.0f} "
                    f"E/I ratio={st['exc_syn']/st['inh_syn'] if st['inh_syn'] else float('inf'):.3f} "
                    f"edges E/I={st['exc_edges']}/{st['inh_edges']} "
                    f"targets={st['n_targets']} fanout mean/max={st['mean_fanout']:.1f}/{st['max_fanout']}"
                )
                print(f"       top types: {st['target_types_top'][:8]}")

        # path products to P1
        print("\n--- signed anatomical path products to P1 ---")
        for label, pop in ((label_a, a), (label_b, b)):
            pp = path_product_to_p1(pop, MAX_HOP)
            print(label, pp["signed_product_by_depth"], "path_hits", pp["path_edge_hits_by_depth"])
            rows[f"{label}_paths"] = pp

        # recurrence
        print("\n--- recurrence (k=3 out-neighborhood ∩ source) ---")
        for label, pop in ((label_a, a), (label_b, b)):
            print(label, recurrence_stats(pop, 3))

        # direct P1
        print("\n--- direct → P1 ---")
        for label, pop in ((label_a, a), (label_b, b)):
            tot = 0.0
            n = 0
            for i in pop:
                tg, w = outs(int(i))
                for t, ww in zip(tg, w):
                    if int(t) in p1set:
                        tot += float(ww)
                        n += 1
            print(label, "direct edges", n, "syn", tot, "sign", signs[pop].tolist()[:20], "...")

        # hop-1 shared vs unique
        t_a = set(rows[label_a][0]["target_set"].tolist())
        t_b = set(rows[label_b][0]["target_set"].tolist())
        print("\n--- hop-1 target overlap ---")
        print("unique_a", len(t_a - t_b), "unique_b", len(t_b - t_a), "shared", len(t_a & t_b))
        print("a-only types", sorted(type_of_set(np.fromiter(t_a - t_b, np.int32)).items(), key=lambda kv: -kv[1])[:15])
        print("b-only types", sorted(type_of_set(np.fromiter(t_b - t_a, np.int32)).items(), key=lambda kv: -kv[1])[:15])
        shared_types = sorted(type_of_set(np.fromiter(t_a & t_b, np.int32)).items(), key=lambda kv: -kv[1])[:15]
        print("shared types", shared_types)

        # who among hop-1 of each touches P1 within 2 more hops (weighted)
        def bridge_score(pop, top_n=15):
            """For each hop-1 target, anatomical product pop→t→...→P1 depth 1..3 from t."""
            acc = []
            for i in pop:
                tg, w = outs(int(i))
                for t, ww in zip(tg, w):
                    t = int(t)
                    if signs[int(i)] == 0:
                        continue
                    # 2-hop from t to P1
                    prod_sum = 0.0
                    tg2, w2 = outs(t)
                    sg = signs[t]
                    if sg == 0:
                        continue
                    for t2, ww2 in zip(tg2, w2):
                        t2 = int(t2)
                        if t2 in p1set:
                            prod_sum += float(ww) * float(ww2) * np.sign(signs[int(i)]) * np.sign(sg)
                        else:
                            tg3, w3 = outs(t2)
                            sg2 = signs[t2]
                            if sg2 == 0:
                                continue
                            for t3, ww3 in zip(tg3, w3):
                                if int(t3) in p1set:
                                    prod_sum += (
                                        float(ww)
                                        * float(ww2)
                                        * float(ww3)
                                        * np.sign(signs[int(i)])
                                        * np.sign(sg)
                                        * np.sign(sg2)
                                    )
                    if prod_sum != 0:
                        acc.append((ntype(t), t, prod_sum, float(ww)))
            # aggregate by type
            by = defaultdict(float)
            for typ, _, p, _ in acc:
                by[typ] += p
            return sorted(by.items(), key=lambda kv: -kv[1])[:top_n], acc

        print(f"\n--- hop1 bridge product to P1 (depth≤3 from source, by type) ---")
        for label, pop in ((label_a, a), (label_b, b)):
            top, _ = bridge_score(pop)
            print(label, top)

        return rows

    rows65 = compare("LgLG6", lg6, "LgLG5", lg5)
    rows78 = compare("LgLG7", lg7, "LgLG8", lg8)

    # save compact json (drop huge target arrays)
    def slim(hops):
        out = []
        for st in hops:
            d = {k: v for k, v in st.items() if k != "target_set"}
            out.append(d)
        return out

    report = {
        "LgLG6_hops": slim(rows65["LgLG6"]),
        "LgLG5_hops": slim(rows65["LgLG5"]),
        "LgLG6_paths": rows65.get("LgLG6_paths"),
        "LgLG5_paths": rows65.get("LgLG5_paths"),
        "LgLG7_hops": slim(rows78["LgLG7"]),
        "LgLG8_hops": slim(rows78["LgLG8"]),
    }
    Path("pathway-decomp-lg6-lg5.json").write_text(json.dumps(report, indent=2, default=str))
    print("\nWrote pathway-decomp-lg6-lg5.json")


if __name__ == "__main__":
    main()
