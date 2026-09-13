"""Build the full annotated MaleCNS v1.0 graph (no connection-weight cutoff).

Standalone copy of the public-data preparation used by the orientation package.
Large source files stay outside git; fetch only when absent.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.feather as feather
from scipy.sparse import csr_matrix, save_npz

SOURCES = [
    (
        "annotations.feather",
        "https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome/body-annotations-male-cns-v1.0-minconf-0.5.feather",
        "2177e246113e4cfbf1e7772ec37c6da1955ff22e8063d0b1f833101f99a9a3b2",
    ),
    (
        "nt.feather",
        "https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome/body-neurotransmitters-male-cns-v1.0.feather",
        "95c9289220663abeb3409f3ad9e5a7f8a53f8093f5139d15502cd08da8879621",
    ),
    (
        "graph.feather",
        "https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome/connectome-weights-male-cns-v1.0-minconf-0.5.feather",
        "e35da783d1c686b2b58b3b87cd6a403ae43bfcfba8bff28e08ef752c1a56afc1",
    ),
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    args = parser.parse_args()
    root = args.data
    root.mkdir(parents=True, exist_ok=True)

    for name, url, digest in SOURCES:
        path = root / name
        if not path.exists():
            temporary = path.with_suffix(".download")
            print("Downloading", url, flush=True)
            urlretrieve(url, temporary)
            temporary.replace(path)
        got = hashlib.sha256(path.read_bytes()).hexdigest()
        if got != digest:
            raise ValueError(f"Source checksum mismatch: {name}")

    rows = feather.read_table(root / "annotations.feather").to_pylist()
    rows = sorted(
        (r for r in rows if r["superclass"] and "tbc" not in r["superclass"]),
        key=lambda r: r["bodyId"],
    )
    ids = np.array([r["bodyId"] for r in rows], dtype=np.int64)
    nttable = feather.read_table(root / "nt.feather", columns=["body", "consensus_nt"])
    nttable = nttable.filter(pc.is_in(nttable["body"], pa.array(ids)))
    nts = dict(zip(nttable["body"].to_pylist(), nttable["consensus_nt"].to_pylist()))
    signs = np.array(
        [
            {"acetylcholine": 1, "gaba": -1, "glutamate": -1}.get(nts.get(int(i)), 0)
            for i in ids
        ],
        dtype=np.int8,
    )
    groups = {
        "candidate_male": [
            i
            for i, r in enumerate(rows)
            if r["type"] in ["LgLG6", "LgLG7"] and r["entryNerve"] == "ProLN"
        ],
        "candidate_female": [
            i
            for i, r in enumerate(rows)
            if r["type"] in ["LgLG5", "LgLG8"] and r["entryNerve"] == "ProLN"
        ],
        "mAL": [
            i
            for i, r in enumerate(rows)
            if (r["type"] or "").startswith("mAL_m") and nts.get(r["bodyId"]) == "gaba"
        ],
        "P1_related": [
            i
            for i, r in enumerate(rows)
            if (r["type"] or "").startswith("pC1_") and r["fruDsx"] == "coexpress_high"
        ],
        "PPN1": [
            i
            for i, r in enumerate(rows)
            if "Kallman 2015: PPN1" in (r["synonyms"] or "")
        ],
    }
    assert all(groups.values()) and len(set(ids)) == len(ids)
    print("Neurons", len(ids), "groups", {k: len(v) for k, v in groups.items()}, flush=True)

    source = pa.memory_map(str(root / "graph.feather"))
    reader = pa.ipc.open_file(source)
    pre, post, weights = [], [], []
    for b in range(reader.num_record_batches):
        t = pa.Table.from_batches([reader.get_batch(b)])
        t = t.filter(
            pc.and_(
                pc.is_in(t["body_pre"], pa.array(ids)),
                pc.is_in(t["body_post"], pa.array(ids)),
            )
        )
        pre.append(np.searchsorted(ids, t["body_pre"].to_numpy()).astype(np.int32))
        post.append(np.searchsorted(ids, t["body_post"].to_numpy()).astype(np.int32))
        weights.append(t["weight"].to_numpy().astype(np.float32))
    pre, post, weights = map(np.concatenate, (pre, post, weights))
    assert np.all(weights > 0)
    graph = csr_matrix((weights, (pre, post)), shape=(len(ids), len(ids)))
    graph.sort_indices()
    save_npz(root / "full-graph.npz", graph)
    np.savez(root / "full-nodes.npz", ids=ids, signs=signs)
    meta = {
        "dataset": "male-cns:v1.0",
        "neurons": len(ids),
        "edges": int(graph.nnz),
        "synapses": int(weights.sum(dtype=np.float64)),
        "nt_counts": dict(Counter(nts.get(int(i)) or "unknown" for i in ids)),
        "zero_fast_weight_neurons": int((signs == 0).sum()),
        "groups": groups,
        "group_annotations": {
            k: [
                {
                    x: rows[i][x]
                    for x in [
                        "bodyId",
                        "type",
                        "receptorType",
                        "synonyms",
                        "entryNerve",
                        "fruDsx",
                    ]
                }
                for i in v
            ]
            for k, v in groups.items()
        },
        "sha256": {name: digest for name, _, digest in SOURCES},
    }
    (root / "full-meta.json").write_text(json.dumps(meta, indent=2))
    print(
        {
            k: v
            for k, v in meta.items()
            if k not in ["groups", "group_annotations", "sha256"]
        },
        flush=True,
    )


if __name__ == "__main__":
    main()
