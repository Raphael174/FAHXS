"""
@ author : Raphaël Aubry

CFD validation study — corrected correlations.
JET-A (POSF10325), Gnielinski coil Nu, Colebrook friction.
Artificial Nu/friction errors tuned to match CFD reference.

Usage (VS Code interactive):
    from hps_combustor.configs import preset_CFD_validation as cfg
    combustor = main_solver(coolantProp=cfg.coolant, hotgasProp=cfg.hotgas,
                            combustorProp=cfg.combustor, numericalProp=cfg.numerical,
                            system_requirements=cfg.sysreqs)
"""

from hps_combustor.input_data import (
    coolantProp, hotgasProp, combustorProp, numericalProp, system_requirements
)

coolant = coolantProp(
    T_out = 396,   # K
)

hotgas = hotgasProp(
    fuel         = "POSF10325",
    p0           = 1e5,    # Pa
    mixing_ratio = 4.121,
    mass_flow_g  = 0.09,   # kg/s
    T_inj_LOX    = 110,    # K
)

combustor = combustorProp(
    inner_diameter      = 0.105,     # m
    exhaust_diameter    = 0.08,      # m
    coil_gap            = 2.18e-3,   # m
    Dh_coil             = 13e-3,     # m
    thickness_coil_wall = 1.5e-3,    # m
    Nusselt_coil        = "Gnielinski",
    friction_coil       = "Colebrook1939",
    material_HX         = "INCO718",
)

numerical = numericalProp(
    L_HX_max                    = 0.57,    # m
    artificial_error_Nu_cold    = -0.85,
    artificial_error_Nu_hot     = -0.51,
    artificial_error_friction_cold = +0.7,
)

sysreqs = system_requirements()
