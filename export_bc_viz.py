"""Export b/c Poisson-stimulation activity under two sign maps for the viewer."""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow.feather as feather

from orientation.experimental import load_ex_network, simulate_ex, type_indices
from orientation.intervention import activate_poisson, none
from orientation.simulate import (
    DURATION_MS,
    apply_sign_overrides,
    load_prepared,
    run_trial,
)

DATA = Path(".experiment-data")
OUT = Path("bc-viz/public/data/trial.json")
SEED = 11
HZ = 80.0
TOP_EDGES = 48


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    net = load_prepared(DATA)
    exnet = load_ex_network(DATA)
    assert np.array_equal(net.ids, exnet.ids)
    ann = feather.read_table(DATA / "annotations.feather").to_pylist()
    by_id = {r["bodyId"]: r for r in ann}
    pos = {
        int(r["id"]): r["position"]
        for r in json.loads((DATA / "full-positions.json").read_text())
    }

    bc = type_indices(exnet, DATA, ["AN09B017b", "AN09B017c"]).astype(np.int32)
    p1 = np.asarray(exnet.groups["P1_readout"], np.int32)
    mal = np.asarray(exnet.groups["mAL"], np.int32)
    lg6 = type_indices(exnet, DATA, ["LgLG6"])
    lg5 = type_indices(exnet, DATA, ["LgLG5"])
    vab3 = np.asarray(exnet.groups["vAB3"], np.int32)

    signs_base = apply_sign_overrides(net)
    signs_exc = signs_base.copy()
    signs_exc[bc] = 1

    interv = activate_poisson(bc, HZ, name="bc_poisson")
    # dedicated RNG stream inside kernel for poisson activation

    print("base WT...", flush=True)
    r_wt = run_trial(net, "no_input", seed=SEED, intervention=none(), signs=signs_base)
    print("base bc stim...", flush=True)
    r_base = run_trial(net, "no_input", seed=SEED, intervention=interv, signs=signs_base)
    print("exc bc stim...", flush=True)
    r_exc = run_trial(net, "no_input", seed=SEED, intervention=interv, signs=signs_exc)

    n = net.n
    has_pos = np.zeros(n, dtype=bool)
    coords = np.zeros((n, 3), dtype=np.float32)
    for i, bid in enumerate(net.ids):
        p = pos.get(int(bid))
        if p is not None:
            has_pos[i] = True
            coords[i] = p
    vis_idx = np.flatnonzero(has_pos)
    c = coords[vis_idx]
    cmin, cmax = c.min(axis=0), c.max(axis=0)
    center = (cmin + cmax) / 2.0
    scale = 200.0 / max(float((cmax - cmin).max()), 1.0)
    coords_n = ((coords - center) * scale).astype(np.float32)
    graph_to_vis = np.full(n, -1, dtype=np.int32)
    graph_to_vis[vis_idx] = np.arange(vis_idx.size, dtype=np.int32)

    bc_set, p1_set, mal_set = set(bc.tolist()), set(p1.tolist()), set(mal.tolist())
    lg6_set, lg5_set, vab3_set = set(lg6.tolist()), set(lg5.tolist()), set(vab3.tolist())

    def role_of(i: int) -> str:
        if i in bc_set:
            return "bc"
        if i in p1_set:
            return "p1"
        if i in mal_set:
            return "mal"
        if i in lg6_set:
            return "lg6"
        if i in lg5_set:
            return "lg5"
        if i in vab3_set:
            return "vab3"
        return "other"

    indptr, indices, data = net.graph.indptr, net.graph.indices, net.graph.data
    edge_acc: dict[int, float] = defaultdict(float)
    for i in bc:
        s, e = indptr[i], indptr[i + 1]
        for t, w in zip(indices[s:e], data[s:e]):
            edge_acc[int(t)] += float(w)
    edges = []
    for t, w in sorted(edge_acc.items(), key=lambda kv: -kv[1])[:TOP_EDGES]:
        v = int(graph_to_vis[t])
        if v >= 0:
            edges.append({"target_vis": v, "weight": float(w)})

    def pack_vis_bins(counts: np.ndarray) -> list:
        out = []
        for b in range(counts.shape[0]):
            act = np.flatnonzero(counts[b])
            pairs = []
            for i in act:
                v = int(graph_to_vis[i])
                if v >= 0:
                    pairs.append([v, int(counts[b, i])])
            out.append(pairs)
        return out

    def p1_spikes(counts: np.ndarray) -> int:
        return int(counts[5:, p1].sum())

    payload = {
        "meta": {
            "seed": SEED,
            "poisson_hz": HZ,
            "duration_ms": float(DURATION_MS),
            "bin_ms": 10.0,
            "n_vis": int(vis_idx.size),
            "conditions": {
                "wt": {"label": "WT · no input", "sign_map": "base"},
                "bc_base": {
                    "label": "b/c Poisson · Glu→−1",
                    "sign_map": "base (b/c inhibitory)",
                },
                "bc_exc": {
                    "label": "b/c Poisson · b/c=+1",
                    "sign_map": "excitatory (vAB3/ACh map)",
                },
            },
            "P1_spikes": {
                "wt": p1_spikes(r_wt.counts),
                "bc_base": p1_spikes(r_base.counts),
                "bc_exc": p1_spikes(r_exc.counts),
            },
            "network_spikes": {
                "wt": int(r_wt.counts.sum()),
                "bc_base": int(r_base.counts.sum()),
                "bc_exc": int(r_exc.counts.sum()),
            },
            "bc_spikes": {
                "wt": int(r_wt.counts[5:, bc].sum()),
                "bc_base": int(r_base.counts[5:, bc].sum()),
                "bc_exc": int(r_exc.counts[5:, bc].sum()),
            },
            "note": "Direct b/c stimulation sensitivity. Preference score not shown.",
        },
        "positions": coords_n[vis_idx].reshape(-1).tolist(),
        "roles": [role_of(int(i)) for i in vis_idx],
        "types": [(by_id.get(int(net.ids[i])) or {}).get("type") or "" for i in vis_idx],
        "bc_vis": [int(graph_to_vis[i]) for i in bc if graph_to_vis[i] >= 0],
        "p1_vis": [int(graph_to_vis[i]) for i in p1 if graph_to_vis[i] >= 0],
        "mal_vis": [int(graph_to_vis[i]) for i in mal if graph_to_vis[i] >= 0],
        "edges": edges,
        "activity": {
            "wt": pack_vis_bins(r_wt.counts),
            "bc_base": pack_vis_bins(r_base.counts),
            "bc_exc": pack_vis_bins(r_exc.counts),
        },
    }
    OUT.write_text(json.dumps(payload, separators=(",", ":")))
    print("wrote", OUT, "bytes", OUT.stat().st_size, flush=True)
    print(json.dumps(payload["meta"], indent=2), flush=True)


if __name__ == "__main__":
    main()
