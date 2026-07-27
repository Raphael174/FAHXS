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
"""

import matplotlib.pyplot as plt
import numpy as np
import CoolProp.CoolProp as CP


_GRID = dict(linewidth=0.8, color='#DDDDDD')
_GRID_MINOR = dict(linestyle='--', linewidth=0.5, color='#EEEEEE')


def _style(ax, x=None):
    ax.set_xlabel(r"$L_{HX}$ [m]")
    if x is not None:
        ax.set_xlim(x[0], x[-1])
    ax.ticklabel_format(useOffset=False)
    ax.grid(which='major', **_GRID)
    ax.grid(which='minor', **_GRID_MINOR)
    ax.minorticks_on()
    ax.legend(fontsize=8)


def _arr(data, key):
    return np.array(data[key])


class HXDashboard:
    def __init__(self, data_master, coolant_name=None):
        self.d = data_master
        self.x = _arr(data_master, "L_HX")
        self.coolant_name = coolant_name or "Coolant"

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
        ax.plot(x, _arr(d, "q_w_rad") / 1e3, color="salmon",  label=r"$q_{w,rad}$")
        ax.set_ylabel(r"Heat flux [kW/m²]")
        _style(ax, x=x)

        ax = axes[1, 0]
        ax.plot(x, d["h_g_conv"], color="orange",          label=r"$h_{g,conv}$")
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

        plt.show()

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

        plt.show()

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

        plt.show()

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

        plt.show()

    # ------------------------------------------------------------------
    def radiation(self):
        """Radiation model outputs."""
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

        plt.show()

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

        plt.show()

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

        plt.show()

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
        ax.plot(x, d["Mach_c"], color="purple",    label=r"$Ma_{He}$")
        ax.plot(x, d["Mach_g"], color="firebrick", label=r"$Ma_g$")
        ax.set_ylabel("Mach [-]")
        _style(ax, x=x)

        # (1,0) Heat flux: total, rad, conv
        ax = axes[1, 0]
        q_w    = _arr(d, "q_w")
        q_rad  = _arr(d, "q_w_rad")
        ax.plot(x, q_w             / 1e3, color="crimson", label=r"$q_w$ total")
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
        ax.plot(x, d["Z"], color="gray", label=r"$Z$")
        ax.set_ylabel(r"Compressibility $Z$ [-]")
        _style(ax, x=x)

        # (3,0) Emissivity / absorptivity
        ax = axes[3, 0]
        ax.plot(x, d["emissivity_g"],   color="firebrick", label=r"$\varepsilon_g$")
        ax.plot(x, d["absorptivity_g"], color="salmon",    label=r"$\alpha_g$")
        ax.set_ylabel("Emissivity / absorptivity [-]")
        _style(ax, x=x)

        # (3,1) Molar fractions CO2 and H2O
        ax = axes[3, 1]
        ax.plot(x, d["X_CO2"], color="gray",      label=r"$X_{CO_2}$")
        ax.plot(x, d["X_H2O"], color="steelblue", label=r"$X_{H_2O}$")
        ax.set_ylabel("Molar fraction [-]")
        _style(ax, x=x)

        # (3,2) Stresses + yield
        ax = axes[3, 2]
        ax.plot(x, np.abs(_arr(d, "stress_inner")) * 1e-6, color="mediumblue", label=r"$|\sigma_{inner}|$")
        ax.plot(x, np.abs(_arr(d, "stress_outer")) * 1e-6, color="darkorange", label=r"$|\sigma_{outer}|$")
        ax.plot(x, _arr(d, "Yield") * 1e-6, color="black", linestyle='--', label="Yield 0.2%")
        ax.set_ylabel("Stress [MPa]")
        _style(ax, x=x)

        plt.show()

    # ------------------------------------------------------------------
    def all(self):
        """Render all thematic figures (boiling()/phase_change() only if quality_c data is present)."""
        self.thermal()
        self.helium()
        self.combustion()
        self.mechanical()
        self.radiation()
        self.boiling()
        self.phase_change()
