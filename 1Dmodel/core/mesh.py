"""Geometry-agnostic 1D flow paths and stream coupling — Stage B of the FV
core rework (docs/solver_design/FV_CORE_REWORK_PLAN.md section 3.1).

The single idea that makes shell-and-tube, shell-and-helical-tube, and both
rocket-nozzle configurations one solver:

    a stream is a 1D path with its OWN arc-length coordinate ``s``, plus a
    monotonic map ``z_of_s`` to the shared axial coordinate ``z`` of the
    assembly.

For a straight tube ``s == z``. For a helical coil, ``s`` runs along the coil
(~1378 nodes for this combustor's real geometry) while ``z`` advances only
one coil pitch per turn. Making ``z_of_s`` explicit DATA computed once by a
geometry builder — rather than a per-step ``dx`` re-derivation inside the
march — is what removes the class of bug the ``HX_config == "shellnHelicalTube"``
guard in ``main_solve.py`` was added to catch (CLAUDE.md, 2026-07-13: the
axial bookkeeping in ``_advance_state()`` silently used a wrong linear-``dx``
approximation for any other config).

``StreamCoupling`` handles the other half: two streams generally do NOT share
a cell partition (one shell cell spans many coil cells), so heat exchanged
between them must be redistributed by a **conservative** overlap operator.
Energy conservation there is exact by construction and is asserted in
tests/test_core_mesh.py — the design doc calls this out as a
"validate before any physics rides on it" item.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FlowPath:
    """One fluid stream's 1D discretization and per-cell geometric measures.

    Conventions, chosen to match the existing solvers so migration is a
    relabeling rather than a reinterpretation:

    - All per-cell areas/perimeters are **per single channel**, not summed
      over parallel channels. ``n_parallel`` carries the multiplicity. (This
      is the opposite of ``transient_core.AxialGrid``, which stored totals —
      the per-channel convention is chosen here because every correlation in
      ``physics/`` is written per-channel, and the existing shell-and-tube
      solver's sharpest edge is exactly the "divide total mass flow by
      N_tubes" bookkeeping. Keeping the per-channel form primary means the
      division happens once, in the builder, instead of at each use site.)
    - ``flow_direction`` is +1 if the fluid flows toward increasing ``z``,
      -1 otherwise. This encodes co- vs counter-flow at the geometry level;
      no solver branch needed.
    - ``s_edges`` is strictly increasing and is the coordinate the fluid
      actually travels along (friction, residence time, entrance length all
      use ``s``). ``z_of_s_edges`` is the projection onto the shared axis
      and must be strictly increasing too (a path that doubles back in ``z``
      is not representable — multi-pass configs need one FlowPath per pass).
    """

    name: str
    s_edges: np.ndarray            # [m] arc length along this stream's own path
    z_of_s_edges: np.ndarray       # [m] shared axial coordinate at each s edge
    A_flow: np.ndarray             # [m^2] per-cell flow area, PER CHANNEL
    Dh: np.ndarray                 # [m] hydraulic diameter
    P_wetted: np.ndarray           # [m] friction perimeter, PER CHANNEL
    P_heated: np.ndarray           # [m] heat-transfer perimeter, PER CHANNEL
    n_parallel: int = 1            # channels/tubes/coil starts in parallel
    geometry_tag: str = "straight_tube"   # closure-registry geometry tag
    orientation_tag: str = "any"          # "vertical" | "horizontal" | "any"
    R_curv: np.ndarray | None = None      # [m] curvature radius; None = straight
    aspect_ratio: np.ndarray | None = None  # h_ch/w_ch for rectangular channels
    inclination: np.ndarray | None = None   # [rad] vs gravity, for buoyancy/HTD
    roughness: float = 1.5e-6      # [m] absolute surface roughness
    flow_direction: int = 1        # +1 toward increasing z, -1 otherwise

    def __post_init__(self) -> None:
        s = _as_1d("s_edges", self.s_edges)
        z = _as_1d("z_of_s_edges", self.z_of_s_edges)
        if s.size < 2:
            raise ValueError("s_edges must contain at least two entries")
        if s.size != z.size:
            raise ValueError("s_edges and z_of_s_edges must have equal length")
        if np.any(np.diff(s) <= 0.0):
            raise ValueError("s_edges must be strictly increasing")
        if np.any(np.diff(z) <= 0.0):
            raise ValueError(
                "z_of_s_edges must be strictly increasing — a path that doubles "
                "back in z needs one FlowPath per pass"
            )
        if self.flow_direction not in (-1, 1):
            raise ValueError("flow_direction must be +1 or -1")
        if int(self.n_parallel) < 1:
            raise ValueError("n_parallel must be >= 1")

        n = s.size - 1
        for field in ("A_flow", "Dh", "P_wetted", "P_heated"):
            arr = _as_cells(field, getattr(self, field), n)
            if np.any(arr <= 0.0):
                raise ValueError(f"{field} must be strictly positive")
            object.__setattr__(self, field, arr)
        for field in ("R_curv", "aspect_ratio", "inclination"):
            val = getattr(self, field)
            if val is not None:
                object.__setattr__(self, field, _as_cells(field, val, n))

        object.__setattr__(self, "s_edges", s)
        object.__setattr__(self, "z_of_s_edges", z)
        object.__setattr__(self, "n_parallel", int(self.n_parallel))

    # -- basic measures -------------------------------------------------
    @property
    def n_cells(self) -> int:
        return self.s_edges.size - 1

    @property
    def ds(self) -> np.ndarray:
        """[m] cell length along the stream's own path (what the fluid sees)."""
        return np.diff(self.s_edges)

    @property
    def dz(self) -> np.ndarray:
        """[m] cell extent projected on the shared axial coordinate."""
        return np.diff(self.z_of_s_edges)

    @property
    def s_centers(self) -> np.ndarray:
        return 0.5 * (self.s_edges[:-1] + self.s_edges[1:])

    @property
    def z_centers(self) -> np.ndarray:
        return 0.5 * (self.z_of_s_edges[:-1] + self.z_of_s_edges[1:])

    @property
    def length_s(self) -> float:
        return float(self.s_edges[-1] - self.s_edges[0])

    @property
    def length_z(self) -> float:
        return float(self.z_of_s_edges[-1] - self.z_of_s_edges[0])

    @property
    def volume_per_channel(self) -> np.ndarray:
        """[m^3] per-cell fluid volume in ONE channel."""
        return self.A_flow * self.ds

    @property
    def volume_total(self) -> np.ndarray:
        """[m^3] per-cell fluid volume summed over all parallel channels."""
        return self.volume_per_channel * self.n_parallel

    @property
    def heated_area_per_channel(self) -> np.ndarray:
        """[m^2] per-cell wetted heat-transfer area in ONE channel."""
        return self.P_heated * self.ds

    @property
    def heated_area_total(self) -> np.ndarray:
        return self.heated_area_per_channel * self.n_parallel

    @property
    def inlet_index(self) -> int:
        return 0 if self.flow_direction == 1 else self.n_cells - 1

    @property
    def outlet_index(self) -> int:
        return self.n_cells - 1 if self.flow_direction == 1 else 0

    def mass_flux(self, mdot_total_kg_s: float) -> np.ndarray:
        """[kg/m^2/s] per-cell mass flux G from the stream's TOTAL mass flow.

        Divides by ``n_parallel`` exactly once, here — this is the
        shell-and-tube "per tube" sharp edge (CLAUDE.md) handled structurally
        instead of at each call site.
        """
        return abs(float(mdot_total_kg_s)) / self.n_parallel / self.A_flow


def _overlap_length(a_lo: float, a_hi: float, b_lo: float, b_hi: float) -> float:
    return max(0.0, min(a_hi, b_hi) - max(a_lo, b_lo))


@dataclass(frozen=True)
class StreamCoupling:
    """Conservative redistribution of an extensive quantity between two
    streams that share an axial axis but not a cell partition.

    ``weights[j, i]`` is the fraction of source cell ``i``'s quantity that
    lands in target cell ``j``, computed from the geometric overlap of their
    ``z`` intervals. For any extensive per-cell array ``q`` (heat rate in W,
    NOT a flux in W/m^2), ``weights @ q`` conserves the sum exactly wherever
    the two streams' ``z`` spans coincide.

    Build with :func:`build_coupling`, which validates the span match.
    """

    weights: np.ndarray  # (n_target, n_source)
    source_name: str
    target_name: str

    def apply(self, q_source: np.ndarray) -> np.ndarray:
        """Redistribute an EXTENSIVE per-cell quantity (e.g. W per cell)."""
        q = np.asarray(q_source, dtype=float)
        if q.shape != (self.weights.shape[1],):
            raise ValueError(
                f"expected source array of shape ({self.weights.shape[1]},), got {q.shape}"
            )
        return self.weights @ q

    @property
    def conservation_defect(self) -> float:
        """Max deviation of each source column's weight sum from 1.

        0 means every source cell's quantity is fully accounted for in the
        target partition (exact conservation). Non-zero means part of the
        source stream's axial extent has no target counterpart.
        """
        return float(np.max(np.abs(self.weights.sum(axis=0) - 1.0)))


def build_coupling(source: FlowPath, target: FlowPath, *, atol: float = 1e-9) -> StreamCoupling:
    """Build the conservative overlap operator mapping ``source`` -> ``target``.

    Both paths must span the same axial interval (within ``atol``); otherwise
    conservation is impossible and this raises rather than silently losing
    energy — the failure mode the design doc explicitly calls out.
    """
    zs = source.z_of_s_edges
    zt = target.z_of_s_edges
    if not (abs(zs[0] - zt[0]) <= atol and abs(zs[-1] - zt[-1]) <= atol):
        raise ValueError(
            f"axial spans differ: {source.name} covers z=[{zs[0]:.6g}, {zs[-1]:.6g}], "
            f"{target.name} covers z=[{zt[0]:.6g}, {zt[-1]:.6g}] — conservative "
            f"coupling requires a common axial interval"
        )

    n_src = source.n_cells
    n_tgt = target.n_cells
    weights = np.zeros((n_tgt, n_src), dtype=float)
    src_dz = source.dz
    for i in range(n_src):
        a_lo, a_hi = zs[i], zs[i + 1]
        for j in range(n_tgt):
            b_lo, b_hi = zt[j], zt[j + 1]
            if b_lo >= a_hi:
                break  # target edges are sorted; no later cell can overlap
            ov = _overlap_length(a_lo, a_hi, b_lo, b_hi)
            if ov > 0.0:
                weights[j, i] = ov / src_dz[i]

    coupling = StreamCoupling(weights=weights, source_name=source.name, target_name=target.name)
    defect = coupling.conservation_defect
    if defect > 1e-10:
        raise ValueError(
            f"overlap operator {source.name}->{target.name} is not conservative "
            f"(max column-sum defect {defect:.3e}); this indicates a gap or "
            f"ordering problem in the axial partitions"
        )
    return coupling


@dataclass(frozen=True)
class HXAssembly:
    """A hot stream, a cold stream, and the conservative maps between them."""

    hot: FlowPath
    cold: FlowPath

    @property
    def hot_to_cold(self) -> StreamCoupling:
        return build_coupling(self.hot, self.cold)

    @property
    def cold_to_hot(self) -> StreamCoupling:
        return build_coupling(self.cold, self.hot)

    @property
    def is_counterflow(self) -> bool:
        return self.hot.flow_direction != self.cold.flow_direction


def _as_1d(name: str, value) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _as_cells(name: str, value, n_cells: int) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        arr = np.full(n_cells, float(arr))
    arr = _as_1d(name, arr)
    if arr.size != n_cells:
        raise ValueError(f"{name} must have shape ({n_cells},), got ({arr.size},)")
    return arr
