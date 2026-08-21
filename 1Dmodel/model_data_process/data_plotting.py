"""
@ author : Raphaël Aubry

HXDashboard — thematic multi-panel figures for VS Code interactive mode.

Usage (in a #%% cell after running the solver):
    db = HXDashboard(combustor.data_master, coolant_name=combustor.coolantProp.coolant)
    db.thermal()        # temperatures and heat transfer
    db.helium()         # coolant flow variables
    db.combustion()     # hot gas flow variables
    db.mechanical()     # stresses and material limits
    db.radiation()      # radiation model outputs
    db.boiling()        # two-phase coolant diagnostics (quality/void/CHF margin,
                         # friction-vs-acceleration dp/dx breakdown); only
                         # meaningful when coolantProp.coolant_model ==
                         # "equilibrium_liquid" (e.g. the Water test config) -
                         # no-op with a message for the single-phase gas march.
    db.phase_change()   # saturation envelope, quality/void, sensible vs latent
                         # heat split - liquid/boiling coolant only.
    db.all()            # all themes in sequence (boiling()/phase_change() only if data present)

Headless / archiving mode:
    db = HXDashboard(data_master, coolant_name=..., save_dir="figures")
    db.all()            # writes one PNG per theme instead of calling plt.show()

    or, as a one-liner (also forces the Agg backend for the duration, then
    restores whatever backend was active):
        save_dashboard(data_master, out_dir, coolant_name=...)

Every panel that depends on a field a given solver does not produce (radiation
outputs, compressibility factor, two-phase fields) is guarded: the panel is
replaced by a short note and the theme still renders. That is what lets the
shell-and-tube adapter (data_plotting_shellntube.py) reuse these exact figures.
"""

import os
from contextlib import contextmanager

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import CoolProp.CoolProp as CP


_GRID = dict(linewidth=0.8, color='#DDDDDD')
_GRID_MINOR = dict(linestyle='--', linewidth=0.5, color='#EEEEEE')


def _style(ax, x=None):
    ax.set_xlabel(r"$L_{HX}$ [m]")
    if x is not None:
        ax.set_xlim(x[0], x[-1])
    try:
        ax.ticklabel_format(useOffset=False)
    except AttributeError:
        pass    # log-scaled axis: LogFormatter has no offset to switch off
    ax.grid(which='major', **_GRID)
    ax.grid(which='minor', **_GRID_MINOR)
    ax.minorticks_on()
    ax.legend(fontsize=8)


def _arr(data, key):
    return np.array(data[key])


def _present(data, key):
    """True when ``key`` holds at least one finite value.

    Distinguishes "this solver does not model that quantity" (absent, empty, or
    all-NaN) from "it does, and here it is" - so one set of panel definitions
    can serve solvers with different physics.
    """
    v = data.get(key)
    if v is None:
        return False
    a = np.asarray(v, dtype=float)
    return a.size > 0 and bool(np.isfinite(a).any())


def _note_absent(ax, msg):
    """Blank a panel whose underlying quantity this run does not model."""
    ax.text(0.5, 0.5, msg, ha='center', va='center', transform=ax.transAxes,
            fontsize=9, color='gray')
    ax.set_xticks([])
    ax.set_yticks([])


@contextmanager
def _agg_backend():
    """Render with Agg, then restore the caller's backend.

    Interactive use of this module must stay unaffected: an archiving call made
    mid-session must not leave the user's own plots headless afterwards.
    """
    previous = matplotlib.get_backend()
    plt.switch_backend("Agg")
    try:
        yield
    finally:
        try:
            plt.switch_backend(previous)
        except Exception:      # pragma: no cover - backend no longer available
            pass


class HXDashboard:
    def __init__(self, data_master, coolant_name=None, save_dir=None, dpi=200):
        self.d = data_master
        self.x = _arr(data_master, "L_HX")
        self.coolant_name = coolant_name or "Coolant"
        self.save_dir = str(save_dir) if save_dir is not None else None
        self.dpi = dpi
        self.written = []
        if self.save_dir is not None:
            os.makedirs(self.save_dir, exist_ok=True)

    # ------------------------------------------------------------------
    def _finish(self, fig, name):
        """Show the figure interactively, or write it into ``save_dir``."""
        if self.save_dir is None:
            plt.show()
            return None
        path = os.path.join(self.save_dir, name)
        fig.savefig(path, dpi=self.dpi, bbox_inches='tight')
        plt.close(fig)
        self.written.append(path)
        return path

    # ------------------------------------------------------------------
    def thermal(self):
        """Temperatures, heat fluxes, heat transfer coefficients, Biot numbers."""
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        fig.suptitle("Thermal overview", fontweight='bold')
        x, d = self.x, self.d

        ax = axes[0, 0]
        ax.plot(x, _arr(d, "T_g") - 273.15,  color="red",            label=r"$T_g$")
        ax.plot(x, _arr(d, "T_wg") - 273.15, color="darkorange",     label=r"$T_{wg}$")
        ax.plot(x, _arr(d, "T_wc") - 273.15, color="mediumblue",     label=r"$T_{wc}$")
        ax.plot(x, _arr(d, "T_c") - 273.15,  color="cornflowerblue", label=r"$T_c$")
        ax.set_ylabel("Temperature [°C]")
        _style(ax, x=x)

        ax = axes[0, 1]
        ax.plot(x, _arr(d, "q_w") / 1e3,     color="crimson", label=r"$q_w$ total")
        if _present(d, "q_w_rad"):
            ax.plot(x, _arr(d, "q_w_rad") / 1e3, color="salmon",  label=r"$q_{w,rad}$")
        ax.set_ylabel(r"Heat flux [kW/m²]")
        _style(ax, x=x)

        ax = axes[1, 0]
        ax.plot(x, d["h_g_conv"], color="orange",          label=r"$h_{g,conv}$")
        if _present(d, "h_g_rad"):
            ax.plot(x, d["h_g_rad"],  color="red",             label=r"$h_{g,rad}$")
        ax.plot(x, d["h_c"],      color="cornflowerblue",  label=r"$h_c$")
        ax.set_ylabel(r"$h$ [W/m²/K]")
        _style(ax, x=x)

        ax = axes[1, 1]
        ax2 = ax.twinx()
        ax.plot(x, d["Biot_g"], color="orange",        label=r"$Bi_g = R_w/R_g$")
        ax.plot(x, d["Biot_c"], color="cornflowerblue", label=r"$Bi_c = R_w/R_c$")
        ax2.plot(x, _arr(d, "UA"), color="green", linestyle='--', label="UA")
        ax.set_ylabel("Biot numbers")
        ax2.set_ylabel("UA [W/K]", color='green')
        ax2.ticklabel_format(useOffset=False)
        lines1, labs1 = ax.get_legend_handles_labels()
        lines2, labs2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labs1 + labs2, fontsize=8)
        _style(ax, x=x)

        self._finish(fig, "thermal.png")

    # ------------------------------------------------------------------
    def helium(self):
        """Coolant flow state along the HX (gas or liquid/boiling)."""
        fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
        fig.suptitle(f"{self.coolant_name} coolant", fontweight='bold')
        x, d = self.x, self.d

        # Sound speed / Mach are real for both gas and liquid modes now (the
        # liquid march uses Wood's equation inside the two-phase dome and the
        # real EOS sound speed outside it - see main_solve.py). Compressibility
        # factor Z is a gas-only concept (no standard liquid/two-phase
        # definition), so its panel is replaced by sound speed - meaningful
        # in both modes rather than all-NaN for liquid coolant.
        panels = [
            (axes[0, 0], _arr(d, "T_c") - 273.15,  "Temperature [°C]",     r"$T_c$",     "cornflowerblue"),
            (axes[0, 1], _arr(d, "p_c") / 1e5,     "Pressure [bar]",       r"$p_c$",     "steelblue"),
            (axes[0, 2], d["rho_c"],                r"Density [kg/m³]",     r"$\rho_c$",  "navy"),
            (axes[1, 0], d["U_c"],                  "Velocity [m/s]",       r"$U_c$",     "teal"),
            (axes[1, 1], d["Mach_c"],               "Mach number",          r"$Ma_c$",    "purple"),
            (axes[1, 2], d["c_c"],                  "Sound speed [m/s]",    r"$c_c$",     "gray"),
        ]
        for ax, y, ylabel, label, color in panels:
            ax.plot(x, y, color=color, label=label)
            ax.set_ylabel(ylabel)
            _style(ax, x=x)

        axes[1, 1].axhline(1.0, color='red', linestyle='--', linewidth=1, label="choking limit")
        axes[1, 1].legend(fontsize=8)

        self._finish(fig, "coolant.png")

    # ------------------------------------------------------------------
    def combustion(self):
        """Hot gas flow state along the HX."""
        fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
        fig.suptitle("Combustion gas", fontweight='bold')
        x, d = self.x, self.d

        panels = [
            (axes[0, 0], _arr(d, "T_g") - 273.15, "Temperature [°C]",   r"$T_g$",     "red"),
            (axes[0, 1], _arr(d, "p_g") / 1e5,    "Pressure [bar]",     r"$p_g$",     "salmon"),
            (axes[0, 2], d["U_g"],                 "Velocity [m/s]",     r"$U_g$",     "darkred"),
            (axes[1, 0], d["Mach_g"],              "Mach number",        r"$Ma_g$",    "firebrick"),
            (axes[1, 1], d["Nu_g"],                "Nusselt (gas side)", r"$Nu_g$",    "orangered"),
            (axes[1, 2], d["cp_g"],                r"$c_p$ [J/kg/K]",   r"$c_{p,g}$", "tomato"),
        ]
        for ax, y, ylabel, label, color in panels:
            ax.plot(x, y, color=color, label=label)
            ax.set_ylabel(ylabel)
            _style(ax, x=x)

        self._finish(fig, "combustion.png")

    # ------------------------------------------------------------------
    def mechanical(self):
        """Coil tube stresses and material limits."""
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        fig.suptitle("Mechanical", fontweight='bold')
        x, d = self.x, self.d

        ax = axes[0, 0]
        ax.plot(x, np.abs(_arr(d, "stress_inner")) * 1e-6, color="mediumblue", label=r"$|\sigma_{inner}|$")
        ax.plot(x, np.abs(_arr(d, "stress_outer")) * 1e-6, color="darkorange", label=r"$|\sigma_{outer}|$")
        ax.plot(x, _arr(d, "Yield") * 1e-6, color="black", linestyle='--', label="Yield 0.2%")
        ax.set_ylabel("Stress [MPa]")
        _style(ax, x=x)

        ax = axes[0, 1]
        ax.plot(x, _arr(d, "stress_thermal_inner") * 1e-6, color="mediumblue", label=r"$\sigma_{th,inner}$")
        ax.plot(x, _arr(d, "stress_thermal_outer") * 1e-6, color="darkorange", label=r"$\sigma_{th,outer}$")
        ax.plot(x, _arr(d, "stress_pressure") * 1e-6,      color="black",      label=r"$\sigma_{pressure}$")
        ax.set_ylabel("Stress [MPa]")
        _style(ax, x=x)

        ax = axes[1, 0]
        yield_arr = _arr(d, "Yield")
        ax.plot(x, np.abs(_arr(d, "stress_inner")) / yield_arr, color="mediumblue", label=r"$|\sigma_{inner}|$ / yield")
        ax.plot(x, np.abs(_arr(d, "stress_outer")) / yield_arr, color="darkorange", label=r"$|\sigma_{outer}|$ / yield")
        ax.axhline(0.8, color='red', linestyle=':', linewidth=1, label="80% limit")
        ax.set_ylabel(r"$\sigma$ / yield")
        _style(ax, x=x)

        ax = axes[1, 1]
        ax.plot(x, _arr(d, "T_wg") - 273.15, color="darkorange", label=r"$T_{wg}$")
        ax.plot(x, _arr(d, "T_wc") - 273.15, color="mediumblue", label=r"$T_{wc}$")
        ax.set_ylabel("Wall temperature [°C]")
        _style(ax, x=x)

        self._finish(fig, "mechanical.png")

    # ------------------------------------------------------------------
    def radiation(self):
        """Radiation model outputs (no-op when the run has no radiation model)."""
        if not _present(self.d, "emissivity_g"):
            print("radiation(): this run has no radiation model outputs "
                  "(emissivity_g/X_CO2/X_H2O absent) - nothing to plot.")
            return
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        fig.suptitle("Radiation", fontweight='bold')
        x, d = self.x, self.d

        ax = axes[0, 0]
        ax.plot(x, d["emissivity_g"],   color="firebrick", label=r"$\varepsilon_g$ (at $T_g$)")
        ax.plot(x, d["absorptivity_g"], color="salmon",    label=r"$\alpha_g$ (at $T_{wg}$)")
        ax.set_ylabel("Gas emissivity / absorptivity")
        _style(ax, x=x)

        ax = axes[0, 1]
        ax.plot(x, _arr(d, "q_w_rad") / 1e3, color="red", label=r"$q_{w,rad}$")
        ax.set_ylabel(r"Radiative flux [kW/m²]")
        _style(ax, x=x)

        ax = axes[1, 0]
        ax.plot(x, d["X_CO2"], color="gray",      label=r"$X_{CO_2}$")
        ax.plot(x, d["X_H2O"], color="steelblue", label=r"$X_{H_2O}$")
        ax.set_ylabel("Molar fractions")
        _style(ax, x=x)

        ax = axes[1, 1]
        h_conv = np.array(d["h_g_conv"], dtype=float)
        h_conv[h_conv == 0] = np.nan
        ax.plot(x, np.array(d["h_g_rad"]) / h_conv, color="red", label=r"$h_{rad}/h_{conv}$")
        ax.set_ylabel(r"$h_{g,rad} / h_{g,conv}$")
        _style(ax, x=x)

        self._finish(fig, "radiation.png")

    # ------------------------------------------------------------------
    def boiling(self):
        """Two-phase coolant diagnostics: quality, void fraction, CHF margin,
        and the friction/acceleration pressure-drop breakdown.

        These fields are only populated when coolantProp.coolant_model ==
        "equilibrium_liquid" (see model_data_process/data_processing.py). For
        the single-phase gas march (e.g. Helium) the arrays are empty; this
        prints a message and returns instead of plotting nothing useful.
        (Enthalpy vs. the saturation envelope is in phase_change() instead.)
        """
        quality = np.asarray(self.d.get("quality_c", []), dtype=float)
        if quality.size == 0:
            print("boiling(): no quality_c/void_c/chf_margin_c data on this run "
                  "- coolantProp.coolant_model was not 'equilibrium_liquid'.")
            return

        fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        fig.suptitle("Boiling coolant diagnostics", fontweight='bold')
        x, d = self.x, self.d

        ax = axes[0, 0]
        ax.plot(x, quality, color="darkorange", label=r"$x$ (quality)")
        ax.axhline(0.0, color='gray', linestyle=':', linewidth=1)
        ax.axhline(1.0, color='red', linestyle=':', linewidth=1, label="x = 1 (dryout)")
        ax.set_ylabel("Vapor quality [-]")
        _style(ax, x=x)

        ax = axes[0, 1]
        ax.plot(x, _arr(d, "void_c"), color="teal", label=r"$\alpha$ (void fraction)")
        ax.set_ylabel("Void fraction [-]")
        _style(ax, x=x)

        ax = axes[1, 0]
        ax.plot(x, _arr(d, "chf_margin_c"), color="crimson", label="CHF margin")
        ax.axhline(1.0, color='red', linestyle='--', linewidth=1, label="CHF limit")
        ax.set_ylabel(r"$q''_{CHF} / q''$ margin [-]")
        _style(ax, x=x)

        ax = axes[1, 1]
        dp_total = _arr(d, "dp_c__dx") / 1e3
        dp_accel = _arr(d, "dp_c__dx_accel") / 1e3
        dp_friction = dp_total - dp_accel
        ax.plot(x, dp_friction, color="steelblue", label="friction")
        ax.plot(x, dp_accel, color="crimson", label="acceleration")
        ax.plot(x, dp_total, color="black", linestyle='--', label="total")
        ax.set_ylabel(r"$dp_c/dx$ [kPa/m]")
        _style(ax, x=x)

        self._finish(fig, "boiling.png")

    # ------------------------------------------------------------------
    def phase_change(self):
        """Phase-change energetics: temperature, enthalpy vs. saturation
        envelope, quality/void, and a sensible-vs-latent heat split.

        The sensible/latent split is a post-hoc diagnostic decomposition of
        the cumulative enthalpy rise, not a separately-tracked solver state:
        each march step's enthalpy change is classified as "latent" if that
        step's quality was inside [0, 1] (two-phase dome) and "sensible"
        otherwise. Only meaningful when coolantProp.coolant_model ==
        "equilibrium_liquid"; no-op with a message otherwise.
        """
        quality = np.asarray(self.d.get("quality_c", []), dtype=float)
        if quality.size == 0:
            print("phase_change(): no quality_c data on this run - "
                  "coolantProp.coolant_model was not 'equilibrium_liquid'.")
            return

        x, d = self.x, self.d
        h = _arr(d, "enthalpy_c")
        p = _arr(d, "p_c")
        Tc = _arr(d, "T_c")

        p_crit = CP.PropsSI("PCRIT", self.coolant_name)
        h_l_sat = np.full_like(p, np.nan)
        h_v_sat = np.full_like(p, np.nan)
        for i in np.where(p < p_crit)[0]:
            h_l_sat[i] = CP.PropsSI("H", "P", p[i], "Q", 0.0, self.coolant_name)
            h_v_sat[i] = CP.PropsSI("H", "P", p[i], "Q", 1.0, self.coolant_name)

        dh = np.diff(h, prepend=h[0])
        in_dome = (quality >= 0.0) & (quality <= 1.0)
        cum_sensible = np.cumsum(np.where(in_dome, 0.0, dh)) / 1e3   # kJ/kg
        cum_latent = np.cumsum(np.where(in_dome, dh, 0.0)) / 1e3     # kJ/kg

        fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        fig.suptitle(f"{self.coolant_name} phase-change energetics", fontweight='bold')

        ax = axes[0, 0]
        ax.plot(x, Tc - 273.15, color="cornflowerblue", label=r"$T_c$")
        ax.set_ylabel("Temperature [°C]")
        _style(ax, x=x)

        ax = axes[0, 1]
        ax.plot(x, h / 1e3, color="mediumblue", label=r"$h_c$")
        ax.plot(x, h_l_sat / 1e3, color="teal", linestyle='--', label=r"$h_{l,sat}(p)$")
        ax.plot(x, h_v_sat / 1e3, color="darkorange", linestyle='--', label=r"$h_{v,sat}(p)$")
        ax.set_ylabel("Enthalpy [kJ/kg]")
        _style(ax, x=x)

        ax = axes[1, 0]
        ax2 = ax.twinx()
        ax.plot(x, quality, color="darkorange", label=r"$x$ (quality)")
        ax2.plot(x, _arr(d, "void_c"), color="teal", linestyle='--', label=r"$\alpha$ (void)")
        ax.axhline(0.0, color='gray', linestyle=':', linewidth=1)
        ax.axhline(1.0, color='red', linestyle=':', linewidth=1)
        ax.set_ylabel("Quality [-]")
        ax2.set_ylabel(r"Void fraction $\alpha$ [-]")
        lines1, labs1 = ax.get_legend_handles_labels()
        lines2, labs2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labs1 + labs2, fontsize=8)
        _style(ax, x=x)

        ax = axes[1, 1]
        ax.plot(x, cum_sensible, color="crimson", label="cumulative sensible heat")
        ax.plot(x, cum_latent, color="purple", label="cumulative latent heat")
        ax.set_ylabel("Cumulative specific heat [kJ/kg]")
        _style(ax, x=x)

        self._finish(fig, "phase_change.png")

    # ------------------------------------------------------------------
    def mega(self):
        """Single-figure overview: 12 panels, all major outputs in one view."""
        fig, axes = plt.subplots(4, 3, figsize=(18, 16), constrained_layout=True)
        fig.suptitle("Full overview", fontweight='bold')
        x, d = self.x, self.d

        def _twin(ax, y_left, y_right, label_left, label_right, color_left, color_right, ylabel_left, ylabel_right):
            ax2 = ax.twinx()
            l1, = ax.plot(x, y_left,  color=color_left,  label=label_left)
            l2, = ax2.plot(x, y_right, color=color_right, linestyle='--', label=label_right)
            ax.set_ylabel(ylabel_left,  color=color_left)
            ax2.set_ylabel(ylabel_right, color=color_right)
            ax2.ticklabel_format(useOffset=False)
            ax.legend(handles=[l1, l2], fontsize=8)
            _style(ax, x=x)

        # (0,0) Temperatures
        ax = axes[0, 0]
        ax.plot(x, _arr(d, "T_g")  - 273.15, color="red",            label=r"$T_g$")
        ax.plot(x, _arr(d, "T_wg") - 273.15, color="darkorange",     label=r"$T_{wg}$")
        ax.plot(x, _arr(d, "T_wc") - 273.15, color="mediumblue",     label=r"$T_{wc}$")
        ax.plot(x, _arr(d, "T_c")  - 273.15, color="cornflowerblue", label=r"$T_c$")
        ax.set_ylabel("Temperature [°C]")
        _style(ax, x=x)

        # (0,1) He pressure
        ax = axes[0, 1]
        ax.plot(x, _arr(d, "p_c") / 1e5, color="steelblue", label=r"$p_c$")
        ax.set_ylabel("Pressure [bar]")
        _style(ax, x=x)

        # (0,2) Mach numbers — He and gas on same axis (both << 1)
        ax = axes[0, 2]
        ax.plot(x, d["Mach_c"], color="purple",    label=r"$Ma_c$")
        ax.plot(x, d["Mach_g"], color="firebrick", label=r"$Ma_g$")
        ax.set_ylabel("Mach [-]")
        _style(ax, x=x)

        # (1,0) Heat flux: total, rad, conv
        ax = axes[1, 0]
        q_w    = _arr(d, "q_w")
        ax.plot(x, q_w             / 1e3, color="crimson", label=r"$q_w$ total")
        if _present(d, "q_w_rad"):
            q_rad  = _arr(d, "q_w_rad")
            ax.plot(x, q_rad           / 1e3, color="salmon",  label=r"$q_{w,rad}$")
            ax.plot(x, (q_w - q_rad)   / 1e3, color="tomato",  label=r"$q_{w,conv}$", linestyle='--')
        ax.set_ylabel(r"Heat flux [kW/m²]")
        _style(ax, x=x)

        # (1,1) Thermal resistances
        ax = axes[1, 1]
        ax.plot(x, d["Res_g"], color="red",            label=r"$R_g$")
        ax.plot(x, d["Res_w"], color="gray",            label=r"$R_w$")
        ax.plot(x, d["Res_c"], color="cornflowerblue",  label=r"$R_c$")
        ax.set_ylabel(r"Thermal resistance [K/W]")
        _style(ax, x=x)

        # (1,2) Reynolds — twin axes (Re_c >> Re_g typically)
        _twin(axes[1, 2],
              y_left=d["Re_c"], y_right=d["Re_g"],
              label_left=r"$Re_c$", label_right=r"$Re_g$",
              color_left="cornflowerblue", color_right="red",
              ylabel_left=r"$Re_c$ [-]", ylabel_right=r"$Re_g$ [-]")

        # (2,0) Nusselt — twin axes (Nu_c >> Nu_g)
        _twin(axes[2, 0],
              y_left=d["Nu_c"], y_right=d["Nu_g"],
              label_left=r"$Nu_c$", label_right=r"$Nu_g$",
              color_left="cornflowerblue", color_right="red",
              ylabel_left=r"$Nu_c$ [-]", ylabel_right=r"$Nu_g$ [-]")

        # (2,1) Biot numbers — twin axes (Biot_c >> Biot_g)
        _twin(axes[2, 1],
              y_left=d["Biot_c"], y_right=d["Biot_g"],
              label_left=r"$Bi_c = R_w/R_c$", label_right=r"$Bi_g = R_w/R_g$",
              color_left="cornflowerblue", color_right="orange",
              ylabel_left=r"$Bi_c$ [-]", ylabel_right=r"$Bi_g$ [-]")

        # (2,2) Compressibility
        ax = axes[2, 2]
        if _present(d, "Z"):
            ax.plot(x, d["Z"], color="gray", label=r"$Z$")
            ax.set_ylabel(r"Compressibility $Z$ [-]")
            _style(ax, x=x)
        else:
            _note_absent(ax, "compressibility factor $Z$\nnot tracked by this solver")

        # (3,0) Emissivity / absorptivity
        ax = axes[3, 0]
        if _present(d, "emissivity_g"):
            ax.plot(x, d["emissivity_g"],   color="firebrick", label=r"$\varepsilon_g$")
            ax.plot(x, d["absorptivity_g"], color="salmon",    label=r"$\alpha_g$")
            ax.set_ylabel("Emissivity / absorptivity [-]")
            _style(ax, x=x)
        else:
            _note_absent(ax, "no radiation model on this run")

        # (3,1) Molar fractions CO2 and H2O
        ax = axes[3, 1]
        if _present(d, "X_CO2"):
            ax.plot(x, d["X_CO2"], color="gray",      label=r"$X_{CO_2}$")
            ax.plot(x, d["X_H2O"], color="steelblue", label=r"$X_{H_2O}$")
            ax.set_ylabel("Molar fraction [-]")
            _style(ax, x=x)
        else:
            _note_absent(ax, "gas composition not tracked per node")

        # (3,2) Stresses + yield
        ax = axes[3, 2]
        ax.plot(x, np.abs(_arr(d, "stress_inner")) * 1e-6, color="mediumblue", label=r"$|\sigma_{inner}|$")
        ax.plot(x, np.abs(_arr(d, "stress_outer")) * 1e-6, color="darkorange", label=r"$|\sigma_{outer}|$")
        ax.plot(x, _arr(d, "Yield") * 1e-6, color="black", linestyle='--', label="Yield 0.2%")
        ax.set_ylabel("Stress [MPa]")
        _style(ax, x=x)

        self._finish(fig, "mega.png")

    # ------------------------------------------------------------------
    def all(self):
        """Render every thematic figure, plus the single-figure ``mega()`` overview.

        Themes whose data this run does not carry (radiation, two-phase) skip
        themselves with a printed note. Returns the list of files written when
        in ``save_dir`` mode, and an empty list interactively.
        """
        self.thermal()
        self.helium()
        self.combustion()
        self.mechanical()
        self.radiation()
        self.boiling()
        self.phase_change()
        self.mega()
        return list(self.written)


def save_dashboard(data_master, out_dir, coolant_name=None, dpi=200):
    """Write the full dashboard to ``out_dir`` as PNGs; return the written paths.

    Forces the Agg backend for the duration so an automated run can never block
    on a GUI window, and restores the previous backend on the way out.
    """
    with _agg_backend():
        db = HXDashboard(data_master, coolant_name=coolant_name,
                         save_dir=out_dir, dpi=dpi)
        return db.all()
