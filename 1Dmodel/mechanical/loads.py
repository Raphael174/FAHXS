""" 
script to compute basic loads on the high pressure high/low temperature helium pipe

Formulas from Roark’s Formulas for Stress and Strain, 7th edition
"""

def stress_pressure_tube (P, thickness_pipe, Dh_pipe):
    """ 
    Roark’s Formulas for Stress and Strain, 7th edition, p562
    returns pressure stress in Pa through pipe thickness
    """
    return P*Dh_pipe/(2*thickness_pipe)

def stress_thermal_tube (T_inner, T_outer, CTE, E, poisson) :
    """
    Roark’s Formulas for Stress and Strain, 7th edition, p761-762
    Returns thermal stresse in Pa at [inner, outer] surfaces
    """
    term = 0.5*CTE*E/(1-poisson)
    # when cold inside and hot outside, +ive tension inside and -ive compression outside
    return [(T_outer-T_inner)*term,  (T_inner-T_outer)*term]


def stress_external_pressure_tube (P_ext, thickness_pipe, Dh_pipe):
    """
    Hoop (membrane) stress in a thin tube under EXTERNAL pressure — Roark's
    Formulas for Stress and Strain, 7th ed., p562, same formula as
    stress_pressure_tube but the load is compressive (negative). Relevant for
    the shell-and-tube config: tubes see high shell-side pressure (~90 bar He)
    outside and low tube-side pressure (~1-5 bar combustion gas) inside — the
    opposite loading direction from the helical coil (internal pressure).
    """
    return -P_ext*Dh_pipe/(2*thickness_pipe)


def collapse_pressure_thin_tube (E, thickness_pipe, Dh_pipe, poisson):
    """
    Elastic (Euler/Bresse-Bryan) critical external pressure for a long thin
    circular tube with unsupported length >> diameter (no stiffening rings),
    Roark's Formulas for Stress and Strain, 7th ed., Table 15.2, case 1:

        P_cr = 2*E / (1-poisson^2) * (t/D)^3

    Returns P_cr [Pa]. Compare |P_ext| / P_cr as a margin (elastic-buckling
    safety factor); a value approaching 1 indicates the tube wall is at risk
    of external-pressure collapse rather than yielding.
    """
    return 2*E/(1-poisson**2) * (thickness_pipe/Dh_pipe)**3