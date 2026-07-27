"""
@ author : Raphaël Aubry

500 kW design point — Shell-and-helical-coil with Inconel 718 coil.
JET-A (POSF10325), O/F = 4.12, INCO718 coil.

Usage (VS Code interactive):
    from hps_combustor.configs import preset_500kW_INCO718 as cfg
    combustor = main_solver(coolantProp=cfg.coolant, hotgasProp=cfg.hotgas,
                            combustorProp=cfg.combustor, numericalProp=cfg.numerical,
                            system_requirements=cfg.sysreqs)
"""

from hps_combustor.input_data import (
    coolantProp, hotgasProp, combustorProp, numericalProp, system_requirements
)

coolant = coolantProp(
    T_out = 700,   # K
)

hotgas = hotgasProp(
    fuel         = "POSF10325",
    p0           = 1e5,    # Pa
    mixing_ratio = 4.121,
    mass_flow_g  = 0.1,    # kg/s
    T_inj_LOX    = 110,    # K
)

combustor = combustorProp(
    inner_diameter      = 0.098,     # m
    coil_gap            = 0.759e-3,  # m
    Dh_coil             = 13.5e-3,   # m
    gap_shell2coil      = 10e-3,     # m
    material_HX         = "INCO718",
)

numerical = numericalProp(
    L_HX_max = 0.7,   # m
)

sysreqs = system_requirements()
