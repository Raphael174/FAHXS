"""
@ author : Raphaël Aubry
"""

from dataclasses import dataclass

@dataclass
class coolantProp :
    coolant : str = "Helium"
    molar_mass : float = 4.002602 # g/mol, https://coolprop.org/fluid_properties/fluids/Helium.html
    mass_flow_c : float = 96e-3 # kg/s 
    T_in : float = 30+273.15 # K 
    T_out : float = 750# K
    p_in : float = 70e5 # Pa 
    p_out : float = 13e5 # Pa 

@dataclass
class hotgasProp : 
    fuel : str = "diesel-C16H34" # "POSF10325"JET A, "gasoline-E5" SP95E5, "gasoline-E10" SP95E10, "diesel-C16H34"
    fuel_density : float = 750 # kg/m3 #https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/672bd2527be152b1d03b0c1f/original/prescreening-of-liquid-density-and-surface-tension-for-synthetic-aviation-turbine-fuels-by-nuclear-magnetic-resonance-atom-types.pdf
    oxidizer : str = "O2" 
    p0 : float = 1.01325e5 # Pa
    mixing_ratio : float = 2.04 #3.827  # O/F 
    mixing_ratio_st : float = 2.85 # OF_st=3.42 in theory but max temp actually occurs at OF=2.85 for Jet (C11H22)
    mass_flow_g : float = 85e-3 # kg/s
    T_g_init : float = 800 # K #* i.e. ignition temperature
    T_inj_LOX : float = 300 # K

@dataclass
class combustorProp:
    HX_config : str = "shellnHelicalTube"   # "coolingjacket" | "coolingcoil" | "shellntube" | "shellnHelicalTube"
    flow_config : str = "co"           # "counter" | "co" — helium flow direction relative to hot gas
    inner_diameter : float = 136e-3         # [m] combustion chamber inner diameter
    mixing_length : float = 0.05            # [m] injection zone length (no HX)
    wall_thickness_inj : float = 5e-3       # [m] injector-plate thickness
    wall_thickness_cc : float = 2e-3        # [m] chamber wall thickness
    exhaust_diameter : float = 136e-3        # [m] nozzle throat diameter (= 30mm - 6mm)
    # shellnHelicalTube geometry
    N_coils : int = 1                       # number of coil passes (currently single-pass)
    coil_gap : float = 5e-3                 # [m] axial gap between consecutive turns
    Dh_coil : float = 9e-3             # [m] coil inner hydraulic diameter
    thickness_coil_wall : float = 2.4e-3    # [m] coil tube wall thickness (250 MPa 316L @ 90 bar, SF 1.15)
    gap_shell2coil : float = 156e-3/2 - 20e-3 - 72e-3/2 - 10e-3/2 - 1e-3    # [m] radial gap between outer coil face and shell wall
    length_2_coil : float = 1e-3            # [m] axial end-clearance between coil and shell ends
    # Correlation selection
    Nusselt_shell : str = "salimpour2008"    # "ahmed_toroid" | "salimpour2008" | "churchill_bernstein_tightcoil" | "churchill_bernstein"
    Nusselt_coil : str = "mori1967"         # "mori1967" | "Gnielinski"
    Nusselt_correction : float = 0.45       # user tuning factor applied on top of Nu_shell
    friction_coil : str = "CurvedPipeAli2024"  # "CurvedPipeAli2024" | "Colebrook1939"
    # Materials
    material_HX : str = "ST316L"           # coil tube material
    material_CC : str = "ST316L"            # chamber shell material
    # Roughness
    combustor_roughness : float = 1.5e-6      # [m] shell-side surface roughness
    channel_roughness : float = 1.5e-6       # [m] coil inner surface roughness

@dataclass
class ToasterProp : 
    burn_time : float = 100 # s 
    p_tank_kerosene : float = 20e5 # Pa
    percentage_vol_pressurant : float = 0.1 # percentage of tank volume for pressurant gas

@dataclass
class numericalProp:
    circular_wall_ON : bool = True      # circular wall for shellnHelicalTube; planar for cooling jacket
    dx : float = 10e-3                  # [m] nominal spatial step (overridden to π*D_coil/N_arc_steps_per_turn for shellnHelicalTube)
    N_arc_steps_per_turn : int = 70      # sub-steps per coil turn; 1 = one full-turn per march step (default)
    dP_max : float = 5e5               # [Pa] max allowable coolant pressure drop
    L_HX_max : float = 0.545 + 50e-3           # [m] max HX axial length
    mass_HX_max : float = 100          # [kg] max HX mass (optimizer constraint)
    T_c_safety_factor : float = 2.5    # lower T_c limit = T_c_min_coolprop * this factor
    radiation_ON : bool = True         # enable WSGGM thermal radiation (shellnHelicalTube only)
    # Chemistry model: "equilibrium" (HP re-equilibrate), "frozen" (no re-equilibrate), "finite_rate" (TODO)
    chemistry_model : str = "equilibrium"
    equilibrium_dh_gas_ON : bool = True  # legacy alias: False => frozen (use chemistry_model instead)
    # Sensitivity / uncertainty perturbations (additive relative errors, 0 = nominal)
    artificial_error_Nu_cold : float = 0.0
    artificial_error_Nu_hot : float = 0.0
    artificial_error_friction_cold : float = 0.0
    # ---- Debug / sanity-check flags ----
    debug_verbose : bool = False           # per-node diagnostic prints inside march loop
    check_energy_balance : bool = True     # ΣdQ vs Δh_g and Δh_c reported at run end
    check_temperature_ordering : bool = True  # warn if T_c > T_wc or T_wg > T_g at any node
    check_mach_limits : bool = True        # warn if Mach_c or Mach_g exceeds 0.5
    check_stress_limits : bool = True      # warn if stress/yield > 0.8 at any node
    yield_at_hot_wall : bool = True        # True → yield at T_wg (conservative); False → yield at mean wall temp
    check_Re_regime : bool = True          # warn if Re_c < 4000 (turbulent correlation in laminar regime)
    check_Z_deviation : bool = True        # warn if |Z - 1| > Z_tolerance
    Z_tolerance : float = 0.05            # compressibility deviation threshold
    energy_balance_tol : float = 0.02     # relative energy imbalance threshold (2 %)


@dataclass
class system_requirements :
    ambient_pressure : float = 101325 # to compute critical pipe size and produced thrust
    max_thrust : float = 1 # N, max thrust for gas exhausts
    burn_time : float = 110 #s


@dataclass
class CorrelationCoefficients:
    """Empirically-uncertain coefficients in the active heat-transfer and friction correlations.

    All values default to the published literature values. Override for calibration against
    experimental data (dp_c, T_g_out).  Identifiability note: with only dp_c as observable,
    ali_c_hi is cleanly identified.  Adding T_g_out also identifies salimpour_a.
    Most sensitive pairs:  (ali_c_hi → dp_c),  (salimpour_a → T_g_out).
    """
    # ---- Shell-side Nu — Salimpour (2008) ----
    # Nu = a * Re^b * Pr^(1/3) * (pitch/D_o)^c * (T_bulk/T_wall)^n
    # Fitted on Pr = 4-15 (water/glycol); extrapolated to Pr ~ 0.65 for combustion gas.
    salimpour_a : float = 0.317    # prefactor                     [Salimpour 2008, Table 1]
    salimpour_b : float = 0.643    # Re exponent
    salimpour_c : float = -0.215   # (coil_pitch / D_o) exponent
    kays_crawford_n : float = 0.25  # gas T-ratio correction exponent [Kays & Crawford 1993]
    # ---- Coil-side Nu — Mori & Nakayama (1967), low-Pr branch (active for He, Pr ~ 0.67) ----
    # Nu = Pr / (a_lo * (Pr^(2/3) - b_lo)) * Re^(4/5) * (d/2R)^(1/10) * (1 + c_lo*Re*(d/2R)^2)^(1/5)
    mori_a_lo : float = 26.2       # denominator scale factor      [Mori & Nakayama 1967, Eq. 14.53a]
    mori_b_lo : float = 0.074      # Pr offset in denominator
    mori_c_lo : float = 0.098      # secondary Dean-number correction
    # ---- Coil friction — Ali et al. (2024) ----
    # I = (Re * alpha^2)^(1/4),  alpha = d/(2*Rc)
    # f = c_lo * I * alpha^(1/2)  for I <= I_split   (weak curvature)
    # f = c_hi * I^(-4/5) * alpha^(1/2)  for I >  I_split   (strong curvature, active at design)
    ali_c_lo : float = 0.316       # low-I  branch prefactor       [Ali et al. 2024]
    ali_c_hi : float = 0.325       # high-I branch prefactor (active at design: I ~ 3)
    ali_I_split : float = 0.868    # Dean-group branch threshold
    # ---- Radiation ----
    # Le = mbl_factor * (V_gas - V_pipe) / A_pipe
    mbl_factor : float = 3.4       # mean beam length factor        [Hottel 1954; 3.6 arbitrary, 3.4 non-grey correction]
    # Hot-side wall emissivity: 0.85 for oxidised/burnt 316L/Inconel [Touloukian & DeWitt 1970]
    # Log-normal prior recommended (strictly positive, multiplicative uncertainty).
    emissivity_wall : float = 0.85 # hot-wall emissivity (burnt/oxidised steel)

