"""
@ author : Raphaël Aubry

Design variable space for shellnHelicalTube HX optimisation.

Each entry:
    "type"      : "continuous" | "categorical" | "bool"
    "range"     : [lower, upper] for continuous, list of options for categorical
    "activated" : 1 = active design variable, 0 = fixed at default
    "target"    : attribute path in the dataclass hierarchy  (dotted: "dataclass.field")
"""

design_variable_space = {
    # --- Hot gas ---
    "mixing_ratio": {
        "type": "continuous", "range": [1.5, 5.0], "activated": 1,
        "target": "hotgasProp.mixing_ratio",
    },
    "mass_flow_g": {
        "type": "continuous", "range": [0.03, 0.30], "activated": 1,
        "target": "hotgasProp.mass_flow_g",
    },

    # --- Combustor shell ---
    "inner_diameter": {
        "type": "continuous", "range": [60e-3, 200e-3], "activated": 1,
        "target": "combustorProp.inner_diameter",
    },
    "wall_thickness_cc": {
        "type": "continuous", "range": [1e-3, 5e-3], "activated": 0,
        "target": "combustorProp.wall_thickness_cc",
    },

    # --- Helical coil geometry ---
    "Dh_coil": {
        "type": "continuous", "range": [6e-3, 25e-3], "activated": 1,
        "target": "combustorProp.Dh_coil",
    },
    "thickness_coil_wall": {
        "type": "continuous", "range": [0.3e-3, 2e-3], "activated": 1,
        "target": "combustorProp.thickness_coil_wall",
    },
    "coil_gap": {
        "type": "continuous", "range": [1e-3, 10e-3], "activated": 1,
        "target": "combustorProp.coil_gap",
    },
    "gap_shell2coil": {
        "type": "continuous", "range": [2e-3, 20e-3], "activated": 1,
        "target": "combustorProp.gap_shell2coil",
    },

    # --- Material ---
    "material_HX": {
        "type": "categorical", "range": ["ST316L", "INCO718", "CuCrZr"], "activated": 0,
        "target": "combustorProp.material_HX",
    },
}
