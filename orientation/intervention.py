"""Neuron-class intervention primitives for orientation experiments.

Interventions act only on outgoing transmission or class drive of a named
neuron set. They never retune global LIF thresholds, resting potential, or
anatomical weights.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

Kind = Literal["none", "output_silence", "activation", "output_gain"]

DEFAULT_MAX_GAIN = 1.0


def _as_index_array(indices: np.ndarray | list[int]) -> np.ndarray:
    arr = np.asarray(indices, dtype=np.int32).reshape(-1)
    if arr.size and np.any(arr < 0):
        raise ValueError("intervention targets must be non-negative indices")
    return np.unique(arr)


@dataclass(frozen=True)
class Intervention:
    """A single class-level intervention.

    kind:
      - none: identity (WT)
      - output_silence: suppress outgoing synaptic transmission
      - activation: tonic steady-state voltage and/or class Poisson drive
      - output_gain: multiply outgoing weights by gain in [0, max_gain]
    """

    name: str
    kind: Kind
    target: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int32))
    tonic_mV: float = 0.0
    activation_hz: float = 0.0
    gain: float = 1.0
    max_gain: float = DEFAULT_MAX_GAIN

    def __post_init__(self) -> None:
        object.__setattr__(self, "target", _as_index_array(self.target))
        if self.kind not in ("none", "output_silence", "activation", "output_gain"):
            raise ValueError(f"unknown intervention kind: {self.kind}")
        if self.kind == "none" and self.target.size:
            raise ValueError("none intervention must have empty target")
        if self.kind == "output_silence" and not self.target.size:
            raise ValueError("output_silence requires a non-empty target class")
        if self.kind == "output_gain":
            if not self.target.size:
                raise ValueError("output_gain requires a non-empty target class")
            if not (0.0 <= self.gain <= self.max_gain + 1e-12):
                raise ValueError(
                    f"gain {self.gain} outside bounded range [0, {self.max_gain}]"
                )
            object.__setattr__(self, "gain", float(min(self.gain, self.max_gain)))
        if self.kind == "activation":
            if not self.target.size:
                raise ValueError("activation requires a non-empty target class")
            if self.tonic_mV < 0.0:
                raise ValueError("tonic_mV must be non-negative")
            if self.activation_hz < 0.0:
                raise ValueError("activation_hz must be non-negative")
            if self.tonic_mV == 0.0 and self.activation_hz == 0.0:
                raise ValueError("activation needs tonic_mV or activation_hz > 0")
        if self.max_gain <= 0:
            raise ValueError("max_gain must be positive")

    @property
    def output_gain_value(self) -> float:
        """Effective outgoing weight multiplier for target neurons."""
        if self.kind == "none":
            return 1.0
        if self.kind == "output_silence":
            return 0.0
        if self.kind == "output_gain":
            return self.gain
        return 1.0

    @property
    def blocks_output(self) -> bool:
        return self.kind == "output_silence" or (
            self.kind == "output_gain" and self.gain == 0.0
        )

    def packed(self, n: int) -> dict[str, np.ndarray | float]:
        """Pack arrays/scalars for the numba kernel."""
        mask = np.zeros(n, dtype=np.bool_)
        if self.target.size:
            if int(self.target.max()) >= n:
                raise ValueError("intervention target index out of range")
            mask[self.target] = True
        if self.kind == "activation":
            tonic = np.zeros(n, dtype=np.float64)
            if self.tonic_mV:
                tonic[self.target] = self.tonic_mV
            return {
                "output_mask": mask,
                "output_gain": self.output_gain_value,
                "tonic": tonic,
                "activation_hz": float(self.activation_hz),
            }
        return {
            "output_mask": mask,
            "output_gain": self.output_gain_value,
            "tonic": np.zeros(n, dtype=np.float64),
            "activation_hz": 0.0,
        }


def none() -> Intervention:
    return Intervention(name="WT", kind="none")


def output_silence(indices: np.ndarray | list[int], name: str | None = None) -> Intervention:
    target = _as_index_array(indices)
    return Intervention(
        name=name or f"silence_n{target.size}",
        kind="output_silence",
        target=target,
        gain=0.0,
    )


def output_gain(
    indices: np.ndarray | list[int],
    gain: float,
    max_gain: float = DEFAULT_MAX_GAIN,
    name: str | None = None,
) -> Intervention:
    target = _as_index_array(indices)
    return Intervention(
        name=name or f"gain_{gain:g}_n{target.size}",
        kind="output_gain",
        target=target,
        gain=float(gain),
        max_gain=float(max_gain),
    )


def activate_tonic(
    indices: np.ndarray | list[int], tonic_mV: float, name: str | None = None
) -> Intervention:
    target = _as_index_array(indices)
    return Intervention(
        name=name or f"tonic_{tonic_mV:g}_n{target.size}",
        kind="activation",
        target=target,
        tonic_mV=float(tonic_mV),
    )


def activate_poisson(
    indices: np.ndarray | list[int], activation_hz: float, name: str | None = None
) -> Intervention:
    target = _as_index_array(indices)
    return Intervention(
        name=name or f"poisson_{activation_hz:g}_n{target.size}",
        kind="activation",
        target=target,
        activation_hz=float(activation_hz),
    )


def activate_tonic_poisson(
    indices: np.ndarray | list[int],
    tonic_mV: float,
    activation_hz: float,
    name: str | None = None,
) -> Intervention:
    target = _as_index_array(indices)
    return Intervention(
        name=name or f"tonic_poisson_n{target.size}",
        kind="activation",
        target=target,
        tonic_mV=float(tonic_mV),
        activation_hz=float(activation_hz),
    )
