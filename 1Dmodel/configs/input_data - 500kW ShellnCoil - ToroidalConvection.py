"""
@ author : Raphaël Aubry

500 kW design point — Shell-and-helical-coil with Ahmed toroidal convection.
Gasoline-E10, O/F = 3.83, ST316L coil.

Usage (VS Code interactive):
    from hps_combustor.configs import preset_500kW_ShellnCoil as cfg
    combustor = main_solver(coolantProp=cfg.coolant, hotgasProp=cfg.hotgas,
                            combustorProp=cfg.combustor, numericalProp=cfg.numerical,
                            system_requirements=cfg.sysreqs)
"""

from hps_combustor.input_data import (
    coolantProp, hotgasProp, combustorProp, numericalProp, system_requirements
)

coolant = coolantProp(
    T_out   = 650,    # K
    p_out   = 85e5,   # Pa
)

hotgas = hotgasProp(
    fuel          = "gasoline-E10",
    p0            = 1e5,       # Pa
    mixing_ratio  = 3.827,
    mass_flow_g   = 0.1,       # kg/s
    T_inj_LOX     = 110,       # K
)

combustor = combustorProp(
    inner_diameter      = 0.1,       # m
    coil_gap            = 0.759e-3,  # m
    Dh_coil             = 13.5e-3,   # m
    thickness_coil_wall = 1.5e-3,    # m
    gap_shell2coil      = 10e-3,     # m
    material_HX         = "ST316L",
)

numerical = numericalProp(
    L_HX_max = 0.7,   # m
)

sysreqs = system_requirements()
