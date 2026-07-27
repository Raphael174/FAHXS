"""
Hot strength properties of steel type 316 S31600
"""
import numpy as np
from scipy.interpolate import interp1d


# STAINLESS STEEL 316L
"""
High temperature characteristics of stainless steels, 
produced by the American Iron and Steel Institute, A designer’s Handbook Series N°9004
"""
##############################################################################################################################################
##############################################################################################################################################


# °C vs Yield Strength 0.2% Offset Pa
#########################################
R02_316L = \
{"x":np.array([27, 149, 260, 371, 482, 593, 704, 816]),
"y" :np.array([290, 201, 172, 159, 148, 140, 131, 110])*1e6}

linear_interp_YieldStrengthR02_316331600 = interp1d(R02_316L["x"], 
                                                    R02_316L["y"], 
                                                    kind='linear')

def YieldStrengthR02_316331600 (T):
    """    
    param T: °C

    return Yield Strength 0.2% Offset in Pa
    """

    if T < 27:
        return 290e6
    elif T>=27 and T<=816:
        return float(linear_interp_YieldStrengthR02_316331600(T))
    else:
        return 110e6

# °C vs Modulus of elasticity Pa
#########################################
ModulusTension316L = \
{"x":np.array([27, 93, 149, 204, 260, 316, 371, 427, 482, 538, 593, 649, 704, 760, 816]),
"y" :np.array([193, 194, 190, 185, 181, 177, 172, 167, 162, 157, 153, 148, 143, 138, 132])*1e9}

linear_interp_ModulusElasticity_316331600 = interp1d(ModulusTension316L["x"], 
                                                    ModulusTension316L["y"], 
                                                    kind='linear')

def ElasticityModulus_316331600 (T):
    """    
    param T: °C

    return Modulus of elasicity in Pa
    """

    if T < 27:
        return 193e9
    elif T>=27 and T<=816:
        return float(linear_interp_ModulusElasticity_316331600(T))
    else:
        return 132e9
    
# CTE 316L
#
#########################################
CTE_316L = \
{"x":np.array([400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700])-273,
"y" :np.array([1.89, 1.917, 1.944, 1.973, 2.002, 2.031, 2.061, 2.092, 2.123, 2.156, 2.188, 2.222, 2.256, 2.291])*1e-5}

linear_interp_CTE_316331600 = interp1d( CTE_316L["x"], 
                                        CTE_316L["y"], 
                                        kind='linear')

def CTE_316331600 (T):
    """    
    param T: °C

    return CTE in 1/C°
    """

    if T < 400-273:
        return 1.89*1e-5
    elif T>=400-273 and T<=1700-273:
        return float(linear_interp_CTE_316331600(T))
    else:
        return 2.291*1e-5

# Lambda 316L
#
#########################################

ST316L= \
{"x":[23.387,155.484,339.355,537.742,651.935,708.065,750.645,788.871],
 "y":[15.301,17.406,19.925,22.556,23.759,24.173,24.361,24.436]}

ST316L_linear_interp = interp1d(np.array(ST316L["x"]), 
                                  ST316L["y"], 
                                  kind='linear')

def compute_ST316L_conductivity (T):
    T_min = 23.387
    T_max = 788.871

    if T < T_min:
        return 15.301
    elif T>=T_min and T<=T_max:
        return float(ST316L_linear_interp(T))
    else:
        return 24.436


# Cp 316L
#   Specific heat [J/kg-K] vs T [°C].
#   Source: AISI "High temperature characteristics of stainless steels",
#   Designer's Handbook Series N°9004 (same reference as the yield/E tables above).
#   Needed by the transient wall-energy ODE (rho*cp*delta*dTbar/dt); irrelevant to
#   the steady solver. See DESIGN_PLAN_shellntube_transient.md section 4.7.
#########################################
Cp_316L = \
{"x": np.array([20, 100, 200, 300, 400, 500, 600, 700, 800]),
 "y": np.array([500, 500, 515, 530, 550, 565, 580, 595, 610])}

Cp_316L_linear_interp = interp1d(Cp_316L["x"], Cp_316L["y"], kind='linear')

def compute_ST316L_cp (T):
    """param T: °C  ->  specific heat [J/kg-K]"""
    if T < 20:
        return 500.0
    elif T <= 800:
        return float(Cp_316L_linear_interp(T))
    else:
        return 610.0


# INCONEL 718
"""
Provide by Abdulah from Mech team
"""


# °C vs Yield Strength 0.2% Offset Pa
#########################################
R02_INCO718 = \
            {   "x": np.array([-240, -129,  -73,  -18,    0,   20,   93,  204,  316,  427,  538,  649,  760]),  
                "y": np.array([1251, 1148, 1096, 1065, 1045, 1034, 1003,  983,  962,  941,  920,  848,  517])*1e6}   

linear_interp_YieldStrengthR02_INCO718 = interp1d(R02_INCO718["x"], 
                                                    R02_INCO718["y"], 
                                                    kind='linear')

def YieldStrengthR02_INCO718 (T):
    """    
    param T: °C

    return Yield Strength 0.2% Offset in Pa
    """

    if T < -240:
        return 1251e6
    elif T>=-240 and T<=760:
        return float(linear_interp_YieldStrengthR02_INCO718(T))
    else:
        return 517e6

# °C vs Modulus of elasticity Pa
#########################################
ModulusTension_INCO718 = \
                        {"x": np.array([-240, -129,  -73,  -18,    0,   20,   93,  204,  316,  427,  538,  649,  760]),  
                        "y": np.array([214.9, 210.8, 208.8, 206.8, 204.7, 202.7, 200.7, 194.6, 188.5, 182.4, 174.3, 166.2, 154.1])*1e9}

linear_interp_ModulusElasticity_INCO718 = interp1d(ModulusTension_INCO718["x"], 
                                                    ModulusTension_INCO718["y"], 
                                                    kind='linear')

def ElasticityModulus_INCO718 (T):
    """    
    param T: °C

    return Modulus of elasicity in Pa
    """

    if T < -240:
        return 214.9e9
    elif T>=-240 and T<=760:
        return float(linear_interp_ModulusElasticity_INCO718(T))
    else:
        return 154.1e9
    
# CTE 
#
#########################################
CTE_INCO718 = \
            {"x": np.array([-240, -129,  -73,  -18,    0,   20,   93,  204,  316,  427,  538,  649,  760]), 
            "y": np.array([ 9.0, 10.8, 11.5, 11.9, 12.2, 12.3, 12.8, 13.3, 13.9, 14.3, 14.5, 14.9, 15.8])*1e-6}

linear_interp_CTE_INCO718 = interp1d(CTE_INCO718["x"], 
                                    CTE_INCO718["y"], 
                                    kind='linear')

def CTE_INCO718 (T):
    """    
    param T: °C

    return CTE in 1/C°
    """

    if T < -240:
        return 9.0e-6
    elif T>=-240 and T<=760:
        return float(linear_interp_CTE_INCO718(T))
    else:
        return 15.8e-6
    
# Lambda 
#
#########################################

PER718= \
{"x": np.array([-18,    0,   20,   93,  204,  316,  427,  538,  649,  760]), 
       "y": np.array([11.1, 11.3, 11.6, 12.7, 14.2, 15.9, 17.5, 19.0, 20.8, 22.3])}

PER718_linear_interp = interp1d(np.array(PER718["x"]), 
                                  PER718["y"], 
                                  kind='linear')

def compute_PER718_conductivity (T):
    """
    T °C
    """
    if T < -18:
        return 11.1
    elif T>=-18 and T<=760:
        return float(PER718_linear_interp(T))
    else:
        #print("Outside range of Inconel718 conductivity data")
        return 22.3


# Cp INCONEL 718
#   Specific heat [J/kg-K] vs T [°C].
#   Source: Special Metals "INCONEL alloy 718" datasheet (SMC-045), Table of
#   Specific Heat. Needed by the transient wall-energy ODE; irrelevant to the
#   steady solver. See DESIGN_PLAN_shellntube_transient.md section 4.7.
#########################################
Cp_INCO718 = \
{"x": np.array([21,  93,  204, 316, 427, 538, 649, 760, 871]),
 "y": np.array([435, 455, 479, 497, 515, 527, 544, 563, 581])}

Cp_INCO718_linear_interp = interp1d(Cp_INCO718["x"], Cp_INCO718["y"], kind='linear')

def compute_INCO718_cp (T):
    """param T: °C  ->  specific heat [J/kg-K]"""
    if T < 21:
        return 435.0
    elif T <= 871:
        return float(Cp_INCO718_linear_interp(T))
    else:
        return 581.0

"""
ASSEMBLE FUNCTIONS FOR QUICK RETURN
"""

# Poisson's ratio — temperature-independent scalars
# Sources: ASME BPVC (316L), AMS 5662 (Inconel 718), ASM Handbook (CuCrZr)
_POISSON = {
    "ST316L":  0.27,
    "INCO718": 0.284,
    "CuCrZr":  0.30,
}


def init_material_temperature_properties(material):
    """
    Returns temperature-dependent property functions and scalar constants for
    the given material.  All functions expect temperature in °C.

    Returns
    -------
    CTE, E, Yield, Lambda : callables  (T_celsius -> property value in SI)
    density                : float  [kg/m³]  (held constant: metals' ρ varies
                             <2 % up to 800 °C, negligible for the wall thermal
                             mass — a documented approximation, not an oversight)
    poisson                : float  [-]
    Cp                     : callable  (T_celsius -> specific heat [J/kg-K]);
                             used only by the transient wall-energy ODE.

    Note: Cp is appended LAST so existing 6-tuple unpackers keep working if they
    ignore the trailing element; the live callers unpack all seven.
    """
    if material == "ST316L":
        CTE     = CTE_316331600
        E       = ElasticityModulus_316331600
        Yield   = YieldStrengthR02_316331600
        Lambda  = compute_ST316L_conductivity
        density = 7.9e3
        Cp      = compute_ST316L_cp
    elif material == "INCO718":
        CTE     = CTE_INCO718
        E       = ElasticityModulus_INCO718
        Yield   = YieldStrengthR02_INCO718
        Lambda  = compute_PER718_conductivity
        density = 8.2e3
        Cp      = compute_INCO718_cp
    else:
        raise ValueError(f"Unknown material '{material}'. Supported: ST316L, INCO718.")

    poisson = _POISSON.get(material, 0.30)
    return CTE, E, Yield, Lambda, density, poisson, Cp
