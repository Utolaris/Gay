"""Signal autopsy: where does male-cue drive die on the way to P1?"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.sparse import load_npz

from orientation.experimental import load_ex_network, simulate_ex, type_indices
from orientation.simulate import apply_sign_overrides, load_prepared, sensory_event_schedule

DATA = Path(".experiment-data")
SEED = 11


def main() -> None:
    net = load_ex_network(DATA)
    pnet = load_prepared(DATA)
    signs = apply_sign_overrides(pnet)
    n = int(net.n)

    def T(*types: str) -> np.ndarray:
        return np.asarray(type_indices(net, DATA, list(types)), dtype=np.int32)

    p1 = np.asarray(net.groups["P1_readout"], np.int32)
    broad = np.asarray(net.groups["P1_related"], np.int32)
    mal = np.asarray(net.groups["mAL"], np.int32)
    lg6 = T("LgLG6")
    lg7 = T("LgLG7")
    lg5 = T("LgLG5")
    a, b, c, d, e, f, g = (
        T("AN09B017a"),
        T("AN09B017b"),
        T("AN09B017c"),
        T("AN09B017d"),
        T("AN09B017e"),
        T("AN09B017f"),
        T("AN09B017g"),
    )
    bcd = np.unique(np.concatenate([b, c, d]))
    bc = np.unique(np.concatenate([b, c]))
    hub035 = T("AN05B035")
    in11a = T("IN05B011a")
    an03 = T("AN03A008")
    ves022 = T("VES022")
    fla = T("FLA001m")
    sip105 = T("SIP105m")
    sip025 = T("SIP025")
    sip100 = T("SIP100m")

    W = load_npz(DATA / "full-graph.npz").tocsr()
    Wcsc = W.tocsc()
    p1i = np.searchsorted(pnet.ids, [12442, 16719, 17867, 20117, 20803, 23968, 519518, 522419])
    p1_exc_src: list[int] = []
    p1_inh_src: list[int] = []
    for t in p1i:
        s0, s1 = Wcsc.indptr[t], Wcsc.indptr[t + 1]
        for k in range(s0, s1):
            p = int(Wcsc.indices[k])
            if signs[p] > 0:
                p1_exc_src.append(p)
            elif signs[p] < 0:
                p1_inh_src.append(p)
    p1_exc = np.unique(p1_exc_src).astype(np.int32)
    p1_inh = np.unique(p1_inh_src).astype(np.int32)
    print("P1 excitatory sources", len(p1_exc), "inhibitory", len(p1_inh), flush=True)
    print("p1_exc max", p1_exc.max() if p1_exc.size else None, "p1_inh max", p1_inh.max() if p1_inh.size else None, flush=True)
    print("all pop maxima:", flush=True)
    for lab, idx in [
        ("lg6", lg6),
        ("b", b),
        ("c", c),
        ("p1", p1),
        ("broad", broad),
        ("mal", mal),
        ("p1_exc", p1_exc),
        ("p1_inh", p1_inh),
        ("ves", ves022),
        ("fla", fla),
    ]:
        print(f"  {lab}: n={len(idx)} max={idx.max() if idx.size else None}", flush=True)

    def spk(counts: np.ndarray, idx) -> int:
        arr = np.atleast_1d(np.asarray(idx, dtype=np.int64))
        if arr.size == 0:
            return 0
        if int(arr.max()) >= int(counts.shape[1]) or int(arr.min()) < 0:
            raise IndexError(
                f"bad index {int(arr.min())}..{int(arr.max())} n={counts.shape[1]} head={arr[:5]}"
            )
        return int(counts[5:, arr.astype(np.int32)].sum())

    pops = [
        ("LgLG6", lg6),
        ("LgLG7", lg7),
        ("AN09B017b", b),
        ("AN09B017c", c),
        ("AN09B017d", d),
        ("AN09B017e", e),
        ("AN09B017f", f),
        ("AN09B017g", g),
        ("AN05B035", hub035),
        ("IN05B011a", in11a),
        ("AN03A008", an03),
        ("mAL", mal),
        ("VES022", ves022),
        ("SIP100m", sip100),
        ("FLA001m", fla),
        ("SIP105m", sip105),
        ("SIP025", sip025),
        ("P1_exc_src", p1_exc),
        ("P1_inh_src", p1_inh),
        ("P1", p1),
        ("broad", broad),
    ]

    def run_pair(gain: np.ndarray, tonic: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
        ton = np.zeros(n) if tonic is None else tonic
        inputs, events = sensory_event_schedule(pnet, "candidate_male", SEED)
        cm, _, _ = simulate_ex(
            pnet.graph.indptr,
            pnet.graph.indices,
            pnet.graph.data,
            signs,
            inputs,
            events,
            gain,
            ton,
            -70.0,
            False,
            True,
            2.0,
            0.02,
            50.0,
            np.zeros(n, np.bool_),
        )
        inputs_f, events_f = sensory_event_schedule(pnet, "candidate_female", SEED)
        cf, _, _ = simulate_ex(
            pnet.graph.indptr,
            pnet.graph.indices,
            pnet.graph.data,
            signs,
            inputs_f,
            events_f,
            gain,
            ton,
            -70.0,
            False,
            True,
            2.0,
            0.02,
            50.0,
            np.zeros(n, np.bool_),
        )
        return cm, cf

    def pack(name: str, cm: np.ndarray, cf: np.ndarray) -> dict:
        m = {}
        for lab, idx in pops:
            try:
                m[lab] = spk(cm, idx)
            except Exception as exc:
                raise RuntimeError(f"{name} pop={lab}") from exc
        m["net"] = int(cm.sum())
        return {
            "cond": name,
            "M": m,
            "F_P1": spk(cf, p1),
            "F_net": int(cf.sum()),
        }

    rows = []

    g = np.ones(n)
    cm, cf = run_pair(g)
    rows.append(pack("WT", cm, cf))

    g = np.ones(n)
    g[mal] = 0
    cm, cf = run_pair(g)
    rows.append(pack("mAL_sil", cm, cf))

    g = np.ones(n)
    g[bc] = 0
    cm, cf = run_pair(g)
    rows.append(pack("bc_sil", cm, cf))

    g = np.ones(n)
    g[mal] = 0
    g[bc] = 0
    cm, cf = run_pair(g)
    rows.append(pack("bc_mAL", cm, cf))

    g = np.ones(n)
    g[mal] = 0
    g[bcd] = 0
    cm, cf = run_pair(g)
    rows.append(pack("bcd_mAL", cm, cf))

    g = np.ones(n)
    g[mal] = 0
    g[in11a] = 0
    cm, cf = run_pair(g)
    rows.append(pack("mAL_IN11a", cm, cf))

    g = np.ones(n)
    g[mal] = 0
    g[sip100] = 0
    cm, cf = run_pair(g)
    rows.append(pack("mAL_SIP100", cm, cf))

    g = np.ones(n)
    g[mal] = 0
    g[ves022] = 0
    cm, cf = run_pair(g)
    rows.append(pack("mAL_VES022", cm, cf))

    g = np.ones(n)
    g[mal] = 0
    ton = np.zeros(n)
    ton[an03] = 12.0
    cm, cf = run_pair(g, ton)
    rows.append(pack("mAL_a008_t12", cm, cf))

    # all P1-inhibitory sources silenced + mAL (upper bound of brake removal)
    g = np.ones(n)
    g[mal] = 0
    g[p1_inh] = 0
    cm, cf = run_pair(g)
    rows.append(pack("mAL_allP1inh", cm, cf))

    keys = [
        "LgLG6",
        "AN09B017b",
        "AN09B017c",
        "AN09B017e",
        "AN09B017f",
        "AN05B035",
        "AN03A008",
        "mAL",
        "VES022",
        "FLA001m",
        "P1_exc_src",
        "P1_inh_src",
        "P1",
        "net",
    ]
    print("\n=== male-cue signal autopsy seed", SEED, "===")
    hdr = f"{'cond':16s} " + " ".join(f"{k[:9]:>9s}" for k in keys) + f" {'F_P1':>5s}"
    print(hdr)
    for r in rows:
        print(
            f"{r['cond'][:16]:16s} "
            + " ".join(f"{r['M'].get(k, 0):9d}" for k in keys)
            + f" {r['F_P1']:5d}"
        )

    Path("male-signal-autopsy.json").write_text(json.dumps(rows, indent=2))
    print("wrote male-signal-autopsy.json")


if __name__ == "__main__":
    main()
