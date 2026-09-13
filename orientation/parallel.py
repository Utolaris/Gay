"""Trial-level process parallelism without changing single-trial dynamics.

Workers share the prepared graph via process-local mmap (or one load from the
same files). The Numba kernel, seed schedule generation, and intervention
packing are unchanged — only the outer trial loop is parallel.
"""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
from scipy.sparse import csr_matrix

from .intervention import Intervention, none, output_gain, output_silence
from .simulate import (
    DEFAULT_REVERSAL_MV,
    Network,
    response_metrics,
    run_trial,
)

# Process-local cache: data_dir -> Network. Graph arrays may be read-only mmap.
_NETWORK_CACHE: dict[str, Network] = {}
_MMAP_DIR_NAME = "graph_mmap"


def default_worker_count(n_trials: int) -> int:
    """Auto workers: physical CPUs, but cap at 4 (best on 8-core Mac for ~18-trial batches)."""
    cpu = os.cpu_count() or 1
    return max(1, min(cpu, n_trials, 4))


def ensure_graph_mmap(data_dir: Path) -> Path:
    """Materialize uncompressed CSR arrays for read-only mmap sharing.

    Writes once under ``.experiment-data/graph_mmap/``. Subsequent workers
    open the ``.npy`` files with ``mmap_mode='r'`` instead of copying the
    compressed NPZ into every process.
    """
    data_dir = Path(data_dir)
    gdir = data_dir / _MMAP_DIR_NAME
    indptr_p = gdir / "indptr.npy"
    indices_p = gdir / "indices.npy"
    data_p = gdir / "data.npy"
    if indptr_p.exists() and indices_p.exists() and data_p.exists():
        return gdir
    from scipy.sparse import load_npz

    gdir.mkdir(parents=True, exist_ok=True)
    graph = load_npz(data_dir / "full-graph.npz").tocsr()
    np.save(indptr_p, graph.indptr)
    np.save(indices_p, graph.indices)
    np.save(data_p, graph.data)
    return gdir


def load_network_shared(data_dir: Path | str) -> Network:
    """Load Network once per process; prefer mmap'd CSR arrays."""
    key = str(Path(data_dir).resolve())
    cached = _NETWORK_CACHE.get(key)
    if cached is not None:
        return cached

    root = Path(data_dir)
    import json

    meta = json.loads((root / "full-meta.json").read_text())
    nodes = np.load(root / "full-nodes.npz")
    ids = np.asarray(nodes["ids"])
    signs = np.asarray(nodes["signs"], dtype=np.int8)

    gdir = ensure_graph_mmap(root)
    indptr = np.load(gdir / "indptr.npy", mmap_mode="r")
    indices = np.load(gdir / "indices.npy", mmap_mode="r")
    data = np.load(gdir / "data.npy", mmap_mode="r")
    n = len(ids)
    # Build CSR then pin the mmap views so workers share file-backed pages
    # instead of keeping private decompressed copies of full-graph.npz.
    graph = csr_matrix(
        (np.asarray(data), np.asarray(indices), np.asarray(indptr)), shape=(n, n)
    )
    try:
        graph.data = data
        graph.indices = indices
        graph.indptr = indptr
    except Exception:
        pass

    groups = {k: np.asarray(v, dtype=np.int32) for k, v in meta["groups"].items()}
    from .simulate import P1_BODY_IDS, VAB3_BODY_IDS

    groups["vAB3"] = np.searchsorted(ids, VAB3_BODY_IDS).astype(np.int32)
    groups["P1_readout"] = np.searchsorted(ids, P1_BODY_IDS).astype(np.int32)
    if not np.array_equal(ids[groups["P1_readout"]], P1_BODY_IDS):
        raise ValueError("P1 body IDs not found in prepared graph")

    net = Network(graph=graph, signs=signs, ids=ids, groups=groups)
    _NETWORK_CACHE[key] = net
    return net


@dataclass(frozen=True)
class TrialSpec:
    """Picklable description of one trial (no large arrays)."""

    data_dir: str
    input_name: str
    seed: int
    intervention_name: str
    intervention_kind: str  # none | output_silence | output_gain | activation
    target: tuple[int, ...] = ()
    gain: float = 1.0
    tonic_mV: float = 0.0
    activation_hz: float = 0.0
    reversal_mV: float = DEFAULT_REVERSAL_MV
    # Optional label for grouping in results
    label: str = ""

    @staticmethod
    def from_intervention(
        data_dir: Path | str,
        input_name: str,
        seed: int,
        intervention: Intervention,
        reversal_mV: float = DEFAULT_REVERSAL_MV,
        label: str = "",
    ) -> "TrialSpec":
        kind = intervention.kind
        if kind == "none":
            kind = "none"
        return TrialSpec(
            data_dir=str(data_dir),
            input_name=input_name,
            seed=int(seed),
            intervention_name=intervention.name,
            intervention_kind=kind,
            target=tuple(int(x) for x in intervention.target),
            gain=float(intervention.gain),
            tonic_mV=float(intervention.tonic_mV),
            activation_hz=float(intervention.activation_hz),
            reversal_mV=float(reversal_mV),
            label=label,
        )


def build_intervention(spec: TrialSpec) -> Intervention:
    if spec.intervention_kind == "none":
        return none()
    if spec.intervention_kind == "output_silence":
        return output_silence(list(spec.target), name=spec.intervention_name)
    if spec.intervention_kind == "output_gain":
        return output_gain(
            list(spec.target),
            spec.gain,
            name=spec.intervention_name,
        )
    if spec.intervention_kind == "activation":
        from .intervention import Intervention as I

        return I(
            name=spec.intervention_name,
            kind="activation",
            target=np.asarray(spec.target, dtype=np.int32),
            tonic_mV=spec.tonic_mV,
            activation_hz=spec.activation_hz,
        )
    raise ValueError(f"unknown intervention_kind {spec.intervention_kind}")


def _limit_numba_threads() -> None:
    """Avoid oversubscription when process pool already uses multiple cores."""
    try:
        from numba import get_num_threads, set_num_threads

        if get_num_threads() != 1:
            set_num_threads(1)
    except Exception:
        pass


def _run_one_trial(spec: TrialSpec) -> dict[str, Any]:
    """Worker entry: load shared network, run one trial, return metrics."""
    _limit_numba_threads()
    network = load_network_shared(spec.data_dir)
    interv = build_intervention(spec)
    result = run_trial(
        network,
        spec.input_name,
        spec.seed,
        interv,
        reversal_mV=spec.reversal_mV,
    )
    metrics = response_metrics(result, network, interv)
    metrics["label"] = spec.label
    return metrics


def run_trial_specs(
    specs: Sequence[TrialSpec],
    workers: int | None = None,
    progress: bool = True,
) -> list[dict[str, Any]]:
    """Run trial specs serially or in a process pool.

    Results are returned in the same order as ``specs``.
    ``workers=1`` uses the in-process path (same as serial baselines).
    """
    specs = list(specs)
    if not specs:
        return []
    n_workers = default_worker_count(len(specs)) if workers is None else max(1, int(workers))

    if n_workers == 1:
        _limit_numba_threads()
        out = []
        for spec in specs:
            row = _run_one_trial(spec)
            if progress:
                import json

                print(json.dumps(row), flush=True)
            out.append(row)
        return out

    out: list[dict[str, Any] | None] = [None] * len(specs)
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_run_one_trial, spec): i for i, spec in enumerate(specs)}
        for fut in as_completed(futures):
            i = futures[fut]
            row = fut.result()
            out[i] = row
            if progress:
                import json

                print(json.dumps(row), flush=True)
    return [r for r in out if r is not None]


def baseline_trial_specs(
    data_dir: Path | str,
    interventions: dict[str, Intervention],
    inputs: Iterable[str],
    seeds: Iterable[int],
    reversal_mV: float = DEFAULT_REVERSAL_MV,
) -> list[TrialSpec]:
    """Cartesian product: intervention × input × seed, in nested loop order."""
    specs: list[TrialSpec] = []
    for input_name in inputs:
        for seed in seeds:
            for name, interv in interventions.items():
                specs.append(
                    TrialSpec.from_intervention(
                        data_dir,
                        input_name,
                        seed,
                        interv,
                        reversal_mV=reversal_mV,
                        label=name,
                    )
                )
    return specs
