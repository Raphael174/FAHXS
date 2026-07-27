"""
@ author : Raphaël Aubry

Maps an optimizer variable vector x → fresh dataclass instances.
Uses the "target" key in design_variable_space to resolve which dataclass field to set.
Never mutates global defaults.
"""

from hps_combustor.input_data import coolantProp, hotgasProp, combustorProp, numericalProp, system_requirements


def build_input_dataclasses(x, design_variable_space):
    """
    Parameters
    ----------
    x : array-like
        Optimizer variable vector (only active variables, in dict-insertion order).
    design_variable_space : dict
        From optimizer_desVar.py — must have "activated" and "target" keys.

    Returns
    -------
    coolant, hotgas, combustor, numerical, sysreqs : fresh dataclass instances
    """
    active = {name: props for name, props in design_variable_space.items()
              if props.get("activated", 0) == 1 and props["type"] != "categorical"}

    if len(x) != len(active):
        raise ValueError(f"Expected {len(active)} active variables, got {len(x)}")

    coolant  = coolantProp()
    hotgas   = hotgasProp()
    combustor = combustorProp()
    numerical = numericalProp()
    sysreqs  = system_requirements()

    obj_map = {
        "coolantProp":     coolant,
        "hotgasProp":      hotgas,
        "combustorProp":   combustor,
        "numericalProp":   numerical,
        "system_requirements": sysreqs,
    }

    for val, (name, props) in zip(x, active.items()):
        target = props.get("target", "")
        if "." not in target:
            raise ValueError(f"Variable '{name}': 'target' must be 'DataclassName.field', got '{target}'")
        dc_name, field = target.split(".", 1)
        obj = obj_map.get(dc_name)
        if obj is None:
            raise ValueError(f"Variable '{name}': unknown dataclass '{dc_name}'")
        if not hasattr(obj, field):
            raise ValueError(f"Variable '{name}': '{dc_name}' has no field '{field}'")
        setattr(obj, field, val)

    return coolant, hotgas, combustor, numerical, sysreqs
