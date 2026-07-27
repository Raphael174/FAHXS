""" 
@ author : Raphaël Aubry
"""

_SOLVER_DATA_KEYS = [
    # Compressibility
    'Z',
    # Thermal resistances
    'Res_g', 'Res_c', 'Res_w',
    # Heat transfer
    'dQ', 'dh_g', 'dq__dx', 'q_w', 'q_w_rad',
    'UP', 'UA',
    # Wall temperatures
    'T_wg', 'T_wc', 'k_w', 'T_c_check',
    # Coolant flow
    'dU_c__dx', 'dT_c__dx', 'dp_c__dx', 'drho_c__dx',
    'cp_c', 'cv_c', 'gamma_c', 'mu_c', 'k_c', 'Pr_c',
    'Re_c', 'Nu_c', 'h_c',
    'f_c', 'f_fd_c',
    'U_c', 'c_c', 'Mach_c',
    'T_c', 'p_c', 'rho_c',
    'De', 'He',
    # Liquid/boiling coolant fields (only populated when
    # coolantProp.coolant_model == "equilibrium_liquid"; see
    # docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md Phase 2)
    'enthalpy_c', 'dh_c__dx', 'quality_c', 'void_c', 'chf_margin_c',
    # Accelerational (HEM) contribution to dp_c__dx, liquid mode only - see
    # main_solve.py's liquid _advance_state() block. dp_c__dx itself already
    # includes this; this field is a diagnostic split-out of the two terms.
    'dp_c__dx_accel',
    # Hot gas flow
    'dp_g__dx',
    'Mach_g', 'c_g', 'W_g',
    'p_g', 'U_g', 'rho_g', 'T_g',
    'cp_g', 'cv_g', 'gamma_g', 'k_g', 'mu_g',
    'Re_g', 'Re_sh', 'Pr_g', 'Nu_g',
    'h_g_conv', 'h_g_rad',
    'emissivity_g', 'absorptivity_g',
    'X_CO2', 'X_H2O',
    # Geometry / lengths
    'L_HX', 'L_ch',
    # Biot numbers
    'Biot_c', 'Biot_g',
    # Tight-coil correction (only populated when Nusselt_shell == "churchill_bernstein_tightcoil")
    'phi_multiplier',
    # Mechanical
    'CTE', 'Modulus', 'Yield',
    'stress_pressure',
    'stress_thermal_inner', 'stress_thermal_outer',
    'stress_inner', 'stress_outer',
]


def make_solver_data():
    """Return a fresh per-run data dictionary. Call once per solver instance."""
    return {k: [] for k in _SOLVER_DATA_KEYS}
