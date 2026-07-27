"""
FPV (Flamelet/Progress-Variable) manifold for the cooling channel — finite-rate
chemistry via a precomputed table, so the transient march makes zero Cantera calls.

@ author : Raphaël Aubry  (finite-rate chemistry track, WP-C1)

See DESIGN_PLAN_shellntube_transient.md section 4b/C1. The cooling channel is a
single already-mixed stream at fixed mixture fraction Z̄ = 1/(1+O/F), so the
thermochemical state is parameterized by two scalars:

  h  : specific enthalpy (removed) — the non-adiabatic axis (wall heat loss)
  Yc : an UNnormalized recombination progress variable,
       Yc = Y_CO2 + Y_H2O - Y_CO   [kg/kg]
       transported by its own net production rate omega_Yc [1/s * kg/kg].

Physics captured: near the hot inlet recombination is fast → local equilibrium;
as the gas cools hard (large enthalpy removed per unit mass at low hot-gas flow)
omega_Yc collapses and the composition FREEZES between the frozen and equilibrium
bounds — the freeze-out that makes low-flow finite-rate matter. Limits recovered
exactly: omega->0 gives frozen; omega->inf gives Yc=Yc_eq(h) = the C0 equilibrium
manifold.

Runtime (in fluid_pass): march Yc alongside h_removed,
    dYc/dx = omega_Yc(h, Yc) / U_g ,
and interpolate {T, rho, mu, k, cp, yH2O, yCO2} from the table. Table lookup uses a
per-h normalized coordinate c = (Yc - Yc_frozen)/(Yc_eq(h) - Yc_frozen) in [0,1]
so the (h, c) grid is regular; Yc itself is the transported (physical) scalar.

Generation is offline and one-time per (O/F, p) — a set of constant-enthalpy
reactor relaxations from the frozen inlet composition toward equilibrium, one per
enthalpy level, sampled onto the (h, c) grid.
"""
import numpy as np
import hashlib
import json
from pathlib import Path


def _Yc_of_Y(Y, iCO2, iH2O, iCO):
    return float(Y[iCO2] + Y[iH2O] - Y[iCO])


def build_fpv_manifold(gas, Y_inlet, T_inlet, p, species_index,
                       n_h=120, n_c=40, T_floor=340.0, t_relax=5e-2, n_t=400,
                       cache_dir=None):
    """
    Build the (h, c) FPV manifold at fixed elemental composition (single Z̄).

    Parameters
    ----------
    gas : a Cantera Solution (mutated locally; caller should not rely on its state
          after this call).
    Y_inlet : inlet burnt-gas mass fractions (the frozen composition, C=0 reference).
    T_inlet, p : inlet temperature [K] and pressure [Pa].
    species_index : dict-like with 'CO2','H2O','CO' -> index, or a callable name->idx.
    n_h, n_c : manifold grid sizes.
    T_floor : cold end of the enthalpy sweep.
    t_relax, n_t : relaxation time [s] and sub-steps for each constant-h reactor path.

    Returns
    -------
    dict manifold with:
      h_grid [n_h]         : specific enthalpy REMOVED (0 at inlet, increasing) [J/kg]
      c_grid [n_c]         : normalized progress 0..1
      Yc_frozen (scalar)   : Yc of the frozen inlet composition
      Yc_eq [n_h]          : equilibrium Yc at each h level
      T,rho,mu,k,cp,xH2O,xCO2,omega_Yc : each [n_h, n_c]
      (omega_Yc is dYc/dt [1/s], the net production rate of the Yc combo)
    """
    if callable(species_index):
        iCO2, iH2O, iCO = species_index('CO2'), species_index('H2O'), species_index('CO')
    else:
        iCO2, iH2O, iCO = species_index['CO2'], species_index['H2O'], species_index['CO']

    cache_path = _cache_path(
        cache_dir, gas, Y_inlet, T_inlet, p, n_h, n_c, T_floor, t_relax, n_t)
    if cache_path is not None and cache_path.exists():
        return _load_cache(cache_path)

    import cantera as ct

    gas.TPY = T_inlet, p, Y_inlet
    h0 = float(gas.enthalpy_mass)
    Yc_frozen = _Yc_of_Y(np.asarray(Y_inlet, float), iCO2, iH2O, iCO)
    MW = gas.molecular_weights  # kg/kmol

    # enthalpy-removed grid (match the C0 equilibrium manifold span)
    h_removed_grid = np.linspace(0.0, h0 - _h_at_Tfloor(gas, p, Y_inlet, T_floor), n_h)
    c_grid = np.linspace(0.0, 1.0, n_c)

    T2 = np.zeros((n_h, n_c)); rho2 = np.zeros((n_h, n_c)); mu2 = np.zeros((n_h, n_c))
    k2 = np.zeros((n_h, n_c)); cp2 = np.zeros((n_h, n_c)); xh2 = np.zeros((n_h, n_c))
    xc2 = np.zeros((n_h, n_c)); om2 = np.zeros((n_h, n_c)); Yc_eq = np.zeros(n_h)

    for jh, hr in enumerate(h_removed_grid):
        h_j = h0 - hr
        # --- true equilibrium state at this enthalpy (the c=1 anchor) ---
        gas.HPY = h_j, p, Y_inlet
        gas.equilibrate('HP')
        Yc_eq_j = _Yc_of_Y(gas.Y, iCO2, iH2O, iCO)
        Yc_eq[jh] = Yc_eq_j
        eq_state = dict(T=gas.T, rho=gas.density, mu=gas.viscosity,
                        k=gas.thermal_conductivity, cp=gas.cp,
                        xh2o=gas.X[iH2O], xco2=gas.X[iCO2], om=0.0)  # omega=0 at equilibrium
        dYc = Yc_eq_j - Yc_frozen
        if abs(dYc) < 1e-12:
            dYc = 1e-12

        # --- constant-enthalpy relaxation from the frozen inlet composition ---
        gas.HPY = h_j, p, Y_inlet
        r = ct.IdealGasConstPressureReactor(gas, energy="on")
        net = ct.ReactorNet([r]); net.rtol = 1e-8; net.atol = 1e-14
        try:
            ph = r.phase
        except AttributeError:
            ph = r.thermo

        traj = {q: [] for q in ("T", "rho", "mu", "k", "cp", "xh2o", "xco2", "om")}
        c_list = [_sample_c(ph, iCO2, iH2O, iCO, Yc_frozen, dYc)]  # c=0 frozen anchor
        _store(ph, MW, iCO2, iH2O, iCO, traj)
        dt = t_relax / n_t
        t = 0.0
        for _ in range(n_t):
            try:
                net.advance(t + dt)
            except Exception:
                # CVODE stiffness failure (common at cold h where chemistry is frozen):
                # stop sampling the trajectory; the c=1 end is anchored to the exact
                # equilibrium state below, so a truncated trajectory is still bracketed.
                break
            t += dt
            c_list.append(_sample_c(ph, iCO2, iH2O, iCO, Yc_frozen, dYc))
            _store(ph, MW, iCO2, iH2O, iCO, traj)
            if c_list[-1] >= 0.999:
                break

        # ALWAYS append the exact equilibrium anchor at c=1 (fixes the cold-h cvode
        # failure that previously left the c=1 edge as a bad extrapolation).
        c_list.append(1.0)
        for q in traj:
            traj[q].append(eq_state[q])

        c_arr = np.clip(np.array(c_list), 0.0, 1.0)
        order = np.argsort(c_arr)          # monotone c for interpolation
        c_arr = c_arr[order]
        for q in traj:
            traj[q] = np.array(traj[q])[order]
        c_u, idx = np.unique(c_arr, return_index=True)

        def _interp(name):
            return np.interp(c_grid, c_u, np.array(traj[name])[idx])
        T2[jh] = _interp("T"); rho2[jh] = _interp("rho"); mu2[jh] = _interp("mu")
        k2[jh] = _interp("k"); cp2[jh] = _interp("cp"); xh2[jh] = _interp("xh2o")
        xc2[jh] = _interp("xco2"); om2[jh] = np.clip(_interp("om"), 0.0, None)

    manifold = dict(h_grid=h_removed_grid, c_grid=c_grid, Yc_frozen=Yc_frozen, Yc_eq=Yc_eq,
                    T=T2, rho=rho2, mu=mu2, k=k2, cp=cp2, xH2O=xh2, xCO2=xc2, omega_Yc=om2)
    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache_path, **manifold)
    return manifold


def _cache_path(cache_dir, gas, Y_inlet, T_inlet, p, n_h, n_c, T_floor, t_relax, n_t):
    if not cache_dir:
        return None
    source = getattr(gas, "source", None) or gas.name
    payload = {
        "source": str(source),
        "species": list(gas.species_names),
        "T_inlet": round(float(T_inlet), 6),
        "p": round(float(p), 3),
        "Y_inlet": np.round(np.asarray(Y_inlet, dtype=float), 12).tolist(),
        "n_h": int(n_h),
        "n_c": int(n_c),
        "T_floor": round(float(T_floor), 6),
        "t_relax": float(t_relax),
        "n_t": int(n_t),
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()[:20]
    return Path(cache_dir) / f"fpv_{digest}.npz"


def _load_cache(path):
    with np.load(path) as data:
        out = {}
        for key in data.files:
            value = data[key]
            out[key] = value.item() if value.shape == () else value.copy()
        return out


def _sample_c(ph, iCO2, iH2O, iCO, Yc_frozen, dYc):
    return (_Yc_of_Y(ph.Y, iCO2, iH2O, iCO) - Yc_frozen) / dYc


def _store(ph, MW, iCO2, iH2O, iCO, traj):
    wdot = ph.net_production_rates
    om = (MW[iCO2] * wdot[iCO2] + MW[iH2O] * wdot[iH2O] - MW[iCO] * wdot[iCO]) / max(ph.density, 1e-9)
    traj["T"].append(ph.T); traj["rho"].append(ph.density); traj["mu"].append(ph.viscosity)
    traj["k"].append(ph.thermal_conductivity); traj["cp"].append(ph.cp)
    traj["xh2o"].append(ph.X[iH2O]); traj["xco2"].append(ph.X[iCO2]); traj["om"].append(om)


def _h_at_Tfloor(gas, p, Y_inlet, T_floor):
    """Enthalpy at the cold-end temperature (frozen composition) for the sweep span."""
    gas.TPY = T_floor, p, Y_inlet
    return float(gas.enthalpy_mass)


class FPVManifold:
    """Runtime interpolator over a built (h, c) manifold. Transports absolute Yc."""

    def __init__(self, m):
        self.m = m
        self.h = m["h_grid"]; self.c = m["c_grid"]
        self.Yc_frozen = m["Yc_frozen"]

    def _c_of(self, h_removed, Yc):
        Yc_eq = float(np.interp(h_removed, self.h, self.m["Yc_eq"]))
        dYc = Yc_eq - self.Yc_frozen
        if abs(dYc) < 1e-12:
            return 0.0
        return float(np.clip((Yc - self.Yc_frozen) / dYc, 0.0, 1.0))

    def _bilin(self, field, h_removed, c):
        ih = np.interp(h_removed, self.h, np.arange(len(self.h)))
        ic = np.interp(c, self.c, np.arange(len(self.c)))
        i0 = int(np.clip(np.floor(ih), 0, len(self.h) - 2)); wi = ih - i0
        j0 = int(np.clip(np.floor(ic), 0, len(self.c) - 2)); wj = ic - j0
        f = field
        return float((f[i0, j0] * (1 - wi) + f[i0 + 1, j0] * wi) * (1 - wj)
                     + (f[i0, j0 + 1] * (1 - wi) + f[i0 + 1, j0 + 1] * wi) * wj)

    def state(self, h_removed, Yc):
        """Return (T, rho, mu, k, cp, xH2O, xCO2, omega_Yc) at (h_removed, Yc)."""
        c = self._c_of(h_removed, Yc)
        m = self.m
        return (self._bilin(m["T"], h_removed, c), self._bilin(m["rho"], h_removed, c),
                self._bilin(m["mu"], h_removed, c), self._bilin(m["k"], h_removed, c),
                self._bilin(m["cp"], h_removed, c), self._bilin(m["xH2O"], h_removed, c),
                self._bilin(m["xCO2"], h_removed, c), self._bilin(m["omega_Yc"], h_removed, c))

    def Yc_inlet(self):
        """Inlet Yc = equilibrium at h_removed=0 (hot end starts at local equilibrium)."""
        return float(self.m["Yc_eq"][0])
