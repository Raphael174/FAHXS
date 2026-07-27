"""
Unsteady 1-D flamelet in mixture-fraction space (Peters 1984/2000 formulation).

Governing equations (unity Lewis number, low-Mach, non-premixed diffusion
flame in mixture-fraction space Z):

    rho dY_k/dt = rho (chi/2) d^2Y_k/dZ^2  +  wdot_k
    rho cp dT/dt = rho cp (chi/2) d^2T/dZ^2
                   + rho (chi/2) [dcp/dZ + sum_k cp_k dY_k/dZ] dT/dZ
                   - sum_k h_k wdot_k  +  dp/dt

Scalar dissipation profile (Peters counterflow form), parameterized by its
stoichiometric value chi_st (the caller's one true forcing input):

    chi(Z) = chi_st * exp(2*[erfc^-1(2*Z_st)]^2 - 2*[erfc^-1(2*Z)]^2)

Numerics: Strang splitting per step
    CN diffusion half-step -> Cantera (CVODE) chemistry full-step -> CN diffusion half-step

Boundary conditions (Dirichlet, two-feed-stream problem):
    Z=0: T_ox,   Y_ox    (oxidizer / lean stream)
    Z=1: T_fuel, Y_fuel  (fuel / rich stream)

This module is a cleaned, standalone extraction of the numerics validated in
HEATV2's rocket-combustion RIF solver (src/hybrid_rocket/physics/chemistry/rif/
flamelet1d.py). It depends on ONLY numpy, scipy.special, and cantera --
no hybrid_rocket import anywhere.

Deliberately removed relative to the source (see ADAPTATION_GUIDE.md
"WHAT WAS REMOVED AND WHY"): spark()/HP-equilibrium re-ignition, the
quench/extinction predicate (is_burning, min-burning-T floor,
consecutive-quench bookkeeping), the reactive-species hot-band chemistry
gating, and all NN-surrogate-shadow / transition-recorder / telemetry
hooks. This kit targets a STEADY or steadily-changing reacting flow with
no ignition/extinction transient to track.

Usage
-----
    from flamelet_kit.flamelet import Flamelet
    fl = Flamelet("gri30.yaml", n_z=65)
    fl.init_mixing(T_ox=300., Y_ox=Y_air, T_fuel=300., Y_fuel=Y_ch4, p=101325.)
    for _ in range(n_steps):
        fl.step(dt=2e-5, p=101325., T_ox=300., Y_ox=Y_air,
                T_fuel=300., Y_fuel=Y_ch4, chi_st=50.)
    print(fl.T_max, fl.Z_st)
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.special import erfcinv

warnings.filterwarnings("ignore", message=".*Temperature.*outside valid range.*")


# ---------------------------------------------------------------------------
# Module helpers
# ---------------------------------------------------------------------------

def _resolve_mechanism(yaml_key: str) -> str:
    """Resolve a mechanism name/path to a loadable Cantera YAML path.

    Accepts (a) a bare Cantera-bundled mechanism name (e.g. "gri30.yaml" or
    "gri30"), which Cantera resolves itself via its own data-path search, or
    (b) an explicit relative/absolute path to a user-supplied YAML file.
    """
    fname = yaml_key if yaml_key.endswith(".yaml") else f"{yaml_key}.yaml"
    p = Path(yaml_key)
    if p.is_absolute() and p.exists():
        return str(p)
    if p.exists():
        return str(p.resolve())
    # Not a local file -- let Cantera try to resolve it against its own
    # bundled-data search path (this is how "gri30.yaml" resolves with zero
    # extra files, which is what example_run.py / the tests rely on).
    return fname


def _build_z_grid(n_z: int, Z_st: float, beta: float = 3.0) -> np.ndarray:
    """
    Two-sided tanh (Roberts) grid clustered around Z_st.

    Left sub-domain [0, Z_st): n_z//2 points, clustered at the Z_st end.
    Right sub-domain [Z_st, 1]: n_z - n_z//2 points, clustered at the Z_st start.
    Guarantees >=15 nodes inside |Z - Z_st| < 0.05 for n_z=65 (verified in tests/).
    """
    n_left = n_z // 2
    n_right = n_z - n_left

    xi_L = np.linspace(0.0, 1.0, n_left + 1)[:-1]  # exclude endpoint (= Z_st)
    Z_L = Z_st * (1.0 + np.sinh(beta * (xi_L - 1.0)) / np.sinh(beta))

    xi_R = np.linspace(0.0, 1.0, n_right)
    Z_R = Z_st + (1.0 - Z_st) * np.sinh(beta * xi_R) / np.sinh(beta)

    Z = np.concatenate([Z_L, Z_R])
    Z[0] = 0.0
    Z[-1] = 1.0
    return Z


def _chi_profile(Z: np.ndarray, Z_st: float, chi_st: float) -> np.ndarray:
    """Peters scalar dissipation profile chi(Z), clipped away from Z in {0,1}.

    chi_st is an ANCHOR value: chi(Z_st) == chi_st by construction. The
    profile's global maximum sits at Z=0.5 (the counterflow's geometric
    mixing-layer center, erfcinv(1)=0), not at Z_st -- this is the correct
    Peters counterflow shape, not a bug (see tests/test_flamelet_kit.py
    test_chi_profile_positive_and_anchored_at_Z_st for the anchor check)."""
    eps = 1e-12
    Z_clip = np.clip(Z, eps, 1.0 - eps)
    A_st = erfcinv(2.0 * np.clip(Z_st, eps, 1.0 - eps)) ** 2
    A_Z = erfcinv(2.0 * Z_clip) ** 2
    return chi_st * np.exp(2.0 * A_st - 2.0 * A_Z)


def _bilger_Z_st(gas, Y_ox: np.ndarray, Y_fuel: np.ndarray,
                  species_names: list) -> float:
    """Stoichiometric mixture fraction via the Bilger (1990) coupling function,
    using elemental mass fractions of C, H, O. Clipped to [0.01, 0.99]."""
    W = np.array([gas.molecular_weights[gas.species_index(sp)]
                  for sp in species_names])

    def beta(Y):
        b = 0.0
        for k, sp in enumerate(species_names):
            idx = gas.species_index(sp)
            comp = gas.species(idx).composition
            n_C = comp.get("C", 0)
            n_H = comp.get("H", 0)
            n_O = comp.get("O", 0)
            b += Y[k] * (2.0 * n_C + 0.5 * n_H - n_O) / W[k]
        return b

    b_ox = beta(Y_ox)
    b_fuel = beta(Y_fuel)
    if abs(b_fuel - b_ox) < 1e-12:
        return 0.14  # fallback
    Z_st = (0.0 - b_ox) / (b_fuel - b_ox)
    return float(np.clip(Z_st, 0.01, 0.99))


# ---------------------------------------------------------------------------
# Tridiagonal Crank-Nicolson solver (Thomas algorithm, interior nodes only)
# ---------------------------------------------------------------------------

def _cn_step(Z: np.ndarray, phi: np.ndarray, D: np.ndarray,
             dt_half: float, bc_left: float, bc_right: float,
             source: Optional[np.ndarray] = None,
             conv: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Crank-Nicolson half-step for d(phi)/dt = D(Z) d^2(phi)/dZ^2 + conv(Z)*d(phi)/dZ + source.

    `conv`: optional first-derivative ("convection") coefficient a(Z), treated
    FULLY IMPLICIT (frozen at substep start = semi-implicit). Used for the
    cp-gradient flamelet term; its explicit treatment overshot on the
    Z_st-clustered grid in the source solver and drove a limit cycle --
    kept fully-implicit here for the same stability reason.

    Dirichlet BCs: phi[0] = bc_left, phi[-1] = bc_right (held fixed).
    Interior nodes 1..n-2 solved by the Thomas algorithm.
    """
    n = len(phi)
    if n < 3:
        return phi.copy()

    h_l = Z[1:] - Z[:-1]
    h_r = h_l.copy()

    hl = h_l[:-1]
    hr = h_r[1:]
    a2 = 2.0 / ((hl + hr) * hl)
    b2 = -2.0 / (hl * hr)
    c2 = 2.0 / ((hl + hr) * hr)

    D_int = D[1:-1]

    sub = -dt_half * D_int[1:] * a2[1:]
    diag = 1.0 - dt_half * D_int * b2
    sup = -dt_half * D_int[:-1] * c2[:-1]

    phi_int = phi[1:-1]
    rhs = (phi_int
           + dt_half * D_int * (a2 * phi[:-2] + b2 * phi_int + c2 * phi[2:]))

    rhs[0] += dt_half * D_int[0] * a2[0] * bc_left
    rhs[-1] += dt_half * D_int[-1] * c2[-1] * bc_right

    if source is not None:
        rhs += dt_half * source[1:-1]

    if conv is not None:
        conv_int = conv[1:-1]
        a1 = -hr / (hl * (hl + hr))
        b1 = (hr - hl) / (hl * hr)
        c1 = hl / (hr * (hl + hr))
        sub = sub - dt_half * conv_int[1:] * a1[1:]
        diag = diag - dt_half * conv_int * b1
        sup = sup - dt_half * conv_int[:-1] * c1[:-1]
        rhs[0] += dt_half * conv_int[0] * a1[0] * bc_left
        rhs[-1] += dt_half * conv_int[-1] * c1[-1] * bc_right

    n_int = len(diag)
    c_ = np.zeros(n_int)
    d_ = np.zeros(n_int)
    c_[0] = sup[0] / diag[0]
    d_[0] = rhs[0] / diag[0]
    for i in range(1, n_int):
        m = diag[i] - sub[i - 1] * c_[i - 1]
        c_[i] = sup[i] / m if i < n_int - 1 else 0.0
        d_[i] = (rhs[i] - sub[i - 1] * d_[i - 1]) / m
    x = np.zeros(n_int)
    x[-1] = d_[-1]
    for i in range(n_int - 2, -1, -1):
        x[i] = d_[i] - c_[i] * x[i + 1]

    phi_new = phi.copy()
    phi_new[1:-1] = x
    phi_new[0] = bc_left
    phi_new[-1] = bc_right
    return phi_new


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class Flamelet:
    """
    Unsteady 1-D flamelet in mixture-fraction space, for a steady or
    steadily-changing (no ignition/extinction) reacting flow.

    Parameters
    ----------
    mechanism : str
        Cantera mechanism YAML name (e.g. "gri30.yaml", bundled with Cantera)
        or a path to a user-supplied mechanism.
    n_z : int
        Number of Z-grid nodes (default 65).
    diff_mask : array-like of float, optional
        Per-species diffusion multiplier (1.0 = normal unity-Le diffusion,
        0.0 = Le -> infinity, i.e. species frozen in Z). Defaults to all-ones
        (unity Lewis number for every species) -- there is no case-specific
        species lumping in this generic kit; supply your own mask if your
        mechanism has a species that should not diffuse in Z.
    cvode_rtol, cvode_atol : float
        Cantera ReactorNet tolerances for the chemistry substep.
    """

    def __init__(
        self,
        mechanism: str,
        n_z: int = 65,
        diff_mask: Optional[np.ndarray] = None,
        cvode_rtol: float = 1.0e-9,
        cvode_atol: float = 1.0e-15,
    ):
        try:
            import cantera as ct
        except ImportError:
            raise ImportError("Cantera is required for the flamelet solver.")

        self._ct = ct
        mech_path = _resolve_mechanism(mechanism)
        self._gas = ct.Solution(mech_path)

        self.n_z = n_z
        self.n_species = self._gas.n_species
        self.species_names = list(self._gas.species_names)
        self._cvode_rtol = float(cvode_rtol)
        self._cvode_atol = float(cvode_atol)

        if diff_mask is not None:
            self.diff_mask = np.asarray(diff_mask, dtype=float)
        else:
            self.diff_mask = np.ones(self.n_species, dtype=float)

        # State arrays (allocated by init_mixing)
        self.Z: Optional[np.ndarray] = None
        self.T: Optional[np.ndarray] = None
        self.Y: Optional[np.ndarray] = None  # (n_z, n_species)
        self.Z_st: Optional[float] = None

        # Boundary conditions (set by init_mixing / step)
        self.T_ox: float = 300.0
        self.T_fuel: float = 300.0
        self.Y_ox: Optional[np.ndarray] = None
        self.Y_fuel: Optional[np.ndarray] = None

        # Pressure state for the dp/dt compression-heating term (relevant
        # if the heat exchanger's operating pressure ramps; if p is truly
        # constant this term is identically zero).
        self._p_prev: float = 1.0e5
        self._dp_dt_ema: float = 0.0
        self._dp_dt_alpha: float = 0.1

        self.step_count: int = 0
        self.t_flamelet: float = 0.0

        self._reactors = None
        self._nets = None

    # ------------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------------

    def init_mixing(
        self,
        T_ox: float,
        Y_ox: np.ndarray,
        T_fuel: float,
        Y_fuel: np.ndarray,
        p: float,
    ) -> None:
        """Set the linear-mixing (frozen, no-reaction) solution as the initial
        condition, build the Z_st-clustered grid, and allocate the persistent
        per-node Cantera reactors."""
        Y_ox = np.asarray(Y_ox, dtype=float)
        Y_fuel = np.asarray(Y_fuel, dtype=float)

        self.T_ox = float(T_ox)
        self.T_fuel = float(T_fuel)
        self.Y_ox = Y_ox.copy()
        self.Y_fuel = Y_fuel.copy()

        self.Z_st = _bilger_Z_st(self._gas, Y_ox, Y_fuel, self.species_names)
        self.Z = _build_z_grid(self.n_z, self.Z_st)

        self.T = T_ox + (T_fuel - T_ox) * self.Z
        self.Y = (Y_ox[np.newaxis, :] * (1.0 - self.Z[:, np.newaxis])
                  + Y_fuel[np.newaxis, :] * self.Z[:, np.newaxis])
        row_sum = self.Y.sum(axis=1, keepdims=True)
        self.Y = np.where(row_sum > 0, self.Y / np.maximum(row_sum, 1e-30), self.Y)

        self._p_prev = float(p)
        self._dp_dt_ema = 0.0
        self._build_reactors(p)

    def _build_reactors(self, p: float) -> None:
        """Create persistent per-node Cantera reactors. Called by init_mixing;
        call again (or use sync_reactors) if T/Y are set directly from outside."""
        ct = self._ct
        self._reactors = []
        self._nets = []
        for i in range(self.n_z):
            g = ct.Solution(self._gas.source)
            g.TPY = float(np.clip(self.T[i], 250.0, 4500.0)), float(p), \
                self._Y_dict(self.Y[i])
            r = ct.IdealGasConstPressureReactor(g, clone=True, name=f"fk_{i}",
                                                 energy="on")
            net = ct.ReactorNet([r])
            net.rtol = self._cvode_rtol
            net.atol = self._cvode_atol
            self._reactors.append(r)
            self._nets.append(net)

    def sync_reactors(self, p: float) -> None:
        """Re-sync the persistent reactors to the current self.T/self.Y state.

        Public utility (renamed from the source's private `_sync_reactors`,
        which existed there only to support `spark()`). Kept here because it
        is generically useful any time external code seeds/overwrites
        `fl.T`/`fl.Y` directly -- e.g. initializing a heat-exchanger flamelet
        bank member from a neighboring already-converged condition instead of
        from cold mixing (see flamelet_bank.py and example_run.py).
        """
        for i, r in enumerate(self._reactors):
            try:
                ph = r.phase
            except AttributeError:
                ph = r.thermo
            ph.TPY = float(np.clip(self.T[i], 250.0, 4500.0)), float(p), \
                self._Y_dict(self.Y[i])
            r.syncState()

    def _Y_dict(self, Y_row: np.ndarray) -> dict:
        d = {}
        for k, sp in enumerate(self.species_names):
            if Y_row[k] > 0.0:
                d[sp] = float(Y_row[k])
        return d if d else {self.species_names[0]: 1.0}

    # ------------------------------------------------------------------
    # Main time-step
    # ------------------------------------------------------------------

    def step(
        self,
        dt: float,
        p: float,
        T_ox: float,
        Y_ox: np.ndarray,
        T_fuel: float,
        Y_fuel: np.ndarray,
        chi_st: float,
    ) -> None:
        """
        Advance the flamelet by dt via Strang splitting:
            CN diffusion half-step -> Cantera chemistry full-step -> CN diffusion half-step.

        `chi_st` is the caller-supplied scalar dissipation at Z_st -- the
        clean decoupling point between this module and whatever closure the
        host application uses to estimate mixing intensity (see
        ADAPTATION_GUIDE.md for heat-exchanger chi_st estimation routes).
        """
        if self.Z is None:
            raise RuntimeError("Call init_mixing() before step().")

        Y_ox = np.asarray(Y_ox, dtype=float)
        Y_fuel = np.asarray(Y_fuel, dtype=float)
        self.T_ox = float(T_ox)
        self.T_fuel = float(T_fuel)
        self.Y_ox = Y_ox.copy()
        self.Y_fuel = Y_fuel.copy()

        if dt > 0:
            dp_dt_raw = (float(p) - self._p_prev) / dt
            self._dp_dt_ema = (self._dp_dt_alpha * dp_dt_raw
                               + (1.0 - self._dp_dt_alpha) * self._dp_dt_ema)
        self._p_prev = float(p)

        chi = _chi_profile(self.Z, self.Z_st, chi_st)
        D = 0.5 * chi

        self._diffusion_step(dt / 2.0, p, D, Y_ox, Y_fuel)
        self._chemistry_step(dt, p)
        self._diffusion_step(dt / 2.0, p, D, Y_ox, Y_fuel)

        self.step_count += 1
        self.t_flamelet += dt

    # ------------------------------------------------------------------
    # Diffusion substep
    # ------------------------------------------------------------------

    def _diffusion_step(
        self,
        dt_half: float,
        p: float,
        D: np.ndarray,
        Y_ox: np.ndarray,
        Y_fuel: np.ndarray,
    ) -> None:
        """CN diffusion half-step for all species and T."""
        for k in range(self.n_species):
            if self.diff_mask[k] < 1e-10:
                continue  # Le -> infinity: no Z-diffusion for this species
            D_k = D * self.diff_mask[k]
            self.Y[:, k] = _cn_step(
                self.Z, self.Y[:, k], D_k, dt_half,
                bc_left=Y_ox[k], bc_right=Y_fuel[k],
            )

        self.Y = np.maximum(self.Y, 0.0)
        row_sum = self.Y.sum(axis=1)
        mask = row_sum > 1e-20
        self.Y[mask] /= row_sum[mask, np.newaxis]

        cp_arr, cp_k_arr, rho_arr = self._compute_cp_rho_fields(p)

        # cp-gradient flamelet term a(Z)*dT/dZ (Pitsch-Peters), folded into the
        # CN operator as a fully-implicit convection coefficient (see _cn_step
        # docstring for why: explicit treatment overshoots on this grid).
        dcp_dZ = np.gradient(cp_arr, self.Z)
        dY_dZ = np.gradient(self.Y, self.Z, axis=0)
        sum_cpk_dYk_dZ = np.sum(cp_k_arr * dY_dZ, axis=1)
        conv = (D / np.maximum(cp_arr, 1.0)) * (dcp_dZ + sum_cpk_dYk_dZ)

        dp_dt_heat = self._dp_dt_ema / np.maximum(rho_arr * cp_arr, 1.0)

        self.T = _cn_step(
            self.Z, self.T, D, dt_half,
            bc_left=self.T_ox, bc_right=self.T_fuel,
            source=dp_dt_heat, conv=conv,
        )
        self.T = np.clip(self.T, 200.0, 4500.0)

    # ------------------------------------------------------------------
    # Chemistry substep
    # ------------------------------------------------------------------

    def _chemistry_step(self, dt: float, p: float) -> None:
        """Advance every interior node with its persistent Cantera reactor.

        Unlike the source solver, ALL interior nodes are always advanced --
        the reactive-species hot-band gating was a rocket-specific cost
        optimization for a transient ignition/extinction problem and is not
        needed (nor correct) for a steady flamelet with no ignition front to
        track.
        """
        for i in range(1, self.n_z - 1):
            r = self._reactors[i]
            net = self._nets[i]
            try:
                ph = r.phase
            except AttributeError:
                ph = r.thermo
            T_old = float(np.clip(self.T[i], 250.0, 4500.0))
            Y_old = np.array(self.Y[i, :], dtype=float, copy=True)
            ph.TPY = T_old, float(p), self._Y_dict(Y_old)
            r.syncState()
            net.initial_time = 0.0
            try:
                net.advance(float(dt))
                self.T[i] = float(np.clip(ph.T, 250.0, 4500.0))
                self.Y[i, :] = np.maximum(np.array(ph.Y, dtype=float), 0.0)
            except Exception:
                pass  # keep current T/Y; chemistry skipped for this node this step

        self.T[0] = self.T_ox
        self.T[-1] = self.T_fuel
        self.Y[0, :] = self.Y_ox
        self.Y[-1, :] = self.Y_fuel

    # ------------------------------------------------------------------
    # Cantera helper: per-node cp and rho
    # ------------------------------------------------------------------

    def _compute_cp_rho_fields(self, p: float):
        """Return cp_arr (n_z,), cp_k_arr (n_z, n_sp), rho_arr (n_z,)."""
        gas = self._gas
        cp_arr = np.zeros(self.n_z)
        cp_k_arr = np.zeros((self.n_z, self.n_species))
        rho_arr = np.zeros(self.n_z)
        for i in range(self.n_z):
            try:
                gas.TPY = float(np.clip(self.T[i], 250.0, 4500.0)), float(p), \
                    self._Y_dict(self.Y[i])
                cp_arr[i] = gas.cp_mass
                cp_k_arr[i] = gas.partial_molar_cp / gas.molecular_weights
                rho_arr[i] = gas.density
            except Exception:
                cp_arr[i] = 1200.0
                cp_k_arr[i] = 1200.0
                rho_arr[i] = p / (280.0 * float(np.clip(self.T[i], 250.0, 4500.0)))
        return cp_arr, cp_k_arr, rho_arr

    # ------------------------------------------------------------------
    # Read-only diagnostics / accessors
    # ------------------------------------------------------------------

    @property
    def T_max(self) -> float:
        return float(np.max(self.T)) if self.T is not None else 0.0

    def T_at_Z(self, Z_query):
        """Interpolate T onto arbitrary Z value(s) (the map-back-to-physical-
        space pattern used by flamelet_bank.field_at)."""
        return np.interp(Z_query, self.Z, self.T)

    def Y_at_Z(self, Z_query, species_name: str):
        k = self.species_names.index(species_name)
        return np.interp(Z_query, self.Z, self.Y[:, k])

    def n_nodes_near_Z_st(self, half_width: float = 0.05) -> int:
        """Number of Z-nodes in [Z_st - half_width, Z_st + half_width]."""
        if self.Z is None:
            return 0
        return int(np.sum(np.abs(self.Z - self.Z_st) <= half_width))
