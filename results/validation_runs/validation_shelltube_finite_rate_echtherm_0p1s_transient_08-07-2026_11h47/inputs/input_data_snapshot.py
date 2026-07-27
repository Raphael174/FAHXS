"""
@ author : Raphaël Aubry
"""

from dataclasses import dataclass
from typing import Optional

@dataclass
class coolantProp :
    coolant : str = "Helium"
    molar_mass : float = 4.002602 # g/mol, https://coolprop.org/fluid_properties/fluids/Helium.html
    mass_flow_c : float = 150e-3 # kg/s 
    T_in : float = 30+273.15 # K 
    T_out : float = 650
    p_in : float = 82e5 # Pa 
    p_out : float = 13e5 # Pa 

@dataclass
class hotgasProp : 
    fuel : str = "diesel-C16H34" # "POSF10325"JET A, "gasoline-E5" SP95E5, "gasoline-E10" SP95E10, "diesel-C16H34"
    fuel_density : float = 750 # kg/m3 #https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/672bd2527be152b1d03b0c1f/original/prescreening-of-liquid-density-and-surface-tension-for-synthetic-aviation-turbine-fuels-by-nuclear-magnetic-resonance-atom-types.pdf
    oxidizer : str = "O2" 
    p0 : float = 1.01325e5 # Pa
    mixing_ratio : float = 3 #3.827  # O/F 
    mixing_ratio_st : float = 2.85 # OF_st=3.42 in theory but max temp actually occurs at OF=2.85 for Jet (C11H22)
    mass_flow_g : float =  100e-3 #85e-3 # kg/s
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
    coil_gap : float = 4e-3                 # [m] axial gap between consecutive turns
    Dh_coil : float = 13.5e-3             # [m] coil inner hydraulic diameter
    thickness_coil_wall : float = 1e-3*0.85    # [m] coil tube wall thickness (250 MPa 316L @ 90 bar, SF 1.15)
    gap_shell2coil : float = 156e-3/2 - 20e-3 - 72e-3/2 - 10e-3/2 - 1e-3      # [m] radial gap between outer coil face and shell wall
    length_2_coil : float = 1e-3            # [m] axial end-clearance between coil and shell ends
    # Correlation selection
    Nusselt_shell : str = "salimpour2008"    # "ahmed_toroid" | "salimpour2008" | "churchill_bernstein_tightcoil" | "churchill_bernstein"
    Nusselt_coil : str = "mori1967"         # "mori1967" | "Gnielinski"
    Nusselt_correction : float = 0.28       # user tuning factor applied on top of Nu_shell
    friction_coil : str = "CurvedPipeAli2024"  # "CurvedPipeAli2024" | "Colebrook1939"
    # Materials
    material_HX : str = "ST316L"           # coil tube material
    material_CC : str = "ST316L"            # chamber shell material
    # Roughness
    combustor_roughness : float = 1.5e-6      # [m] shell-side surface roughness
    channel_roughness : float = 1.5e-6       # [m] coil inner surface roughness

@dataclass
class shellTubeProp:
    """Baffled shell-and-tube config (hot gas in straight tubes, coolant in shell-side
    baffled cross-flow). Geometry mirrors the EchTherm GEOMETRY screen. See
    DESIGN_PLAN_shellntube_transient.md sections 2-3."""
    # tubes
    N_tubes : int = 235
    D_tube_outer : float = 5e-3
    thickness_tube_wall : float = 0.75e-3
    L_tube : float = 235e-3
    layout : str = "triangular30"        # "triangular30" | "square90" | "rotated45"
    pitch_ratio : float = 1.3
    # shell
    D_shell_inner : float = 110e-3
    shell_thickness : float = 8e-3
    # baffles (single segmental, tubes in windows)
    N_baffles : int = 15
    baffle_cut : float = 0.20            # window opening as fraction of D_shell
    baffle_thickness : float = 3e-3
    baffle_spacing : float = 12e-3       # central baffle spacing shown by EchTherm
    L_front_end : float = 100e-3         # tube-sheet/nozzle end zone, not a baffle spacing
    L_rear_end : float = 10e-3           # tube-sheet/nozzle end zone, not a baffle spacing
    L_inlet_spacing : float = 8e-3       # first baffle spacing from tube inlet
    L_outlet_spacing : float = 8e-3      # last baffle spacing to tube outlet
    tube_sheet_thickness : float = 3e-3
    # diametral clearances + bypass
    clearance_tube_baffle : float = 1e-3
    clearance_baffle_shell : float = 1e-3
    clearance_bundle_shell : float = 0e-3
    N_sealing_strip_pairs : int = 0
    # tube/nozzle surface options from EchTherm geometry page
    inside_tube_choice : str = "grooved"     # "smooth" | "grooved" | "helical_insert" | "intensification_factor" | "power_law"
    outside_tube_choice : str = "smooth"     # "smooth" | "low_finned"
    corrugation_thickness : float = 0.2e-3
    corrugation_pitch : float = 2e-3
    corrugation_angle_deg : float = 30.0
    corrugation_sharp_corner_number : int = 2
    corrugation_sharp_angle_deg : float = 60.0
    tube_side_nozzle_diameter_in : float = 100e-3
    tube_side_nozzle_diameter_out : float = 100e-3
    shell_side_nozzle_diameter_in : float = 21e-3
    shell_side_nozzle_diameter_out : float = 21e-3
    tube_side_nozzle_orientation_in : str = "axial"       # "axial" | "transversal"
    tube_side_nozzle_orientation_out : str = "axial"
    shell_side_nozzle_orientation_in : str = "axial"
    shell_side_nozzle_orientation_out : str = "axial"
    # fluid allocation + correlation selection + materials
    tube_side_fluid : str = "hotgas"     # "hotgas" | "coolant" (hotgas implemented first)
    Nusselt_tube : str = "gnielinski_blended"
    Nusselt_shell_baffled : str = "bell_delaware"
    material_tube : str = "INCO718"
    material_shell : str = "INCO718"
    tube_roughness : float = 1.5e-6
    inlet_mode : str = "combustor"       # "combustor" (fed by Cantera combustion) | "prescribed"


@dataclass
class ToasterProp :
    burn_time : float = 100 # s 
    p_tank_kerosene : float = 20e5 # Pa
    percentage_vol_pressurant : float = 0.1 # percentage of tank volume for pressurant gas

@dataclass
class numericalProp:
    circular_wall_ON : bool = True      # circular wall for shellnHelicalTube; planar for cooling jacket
    dx : float = 10e-3                  # [m] nominal spatial step (overridden to π*D_coil/N_arc_steps_per_turn for shellnHelicalTube)
    N_arc_steps_per_turn : int = 50      # sub-steps per coil turn; 1 = one full-turn per march step (default)
    dP_max : float = 5e5               # [Pa] max allowable coolant pressure drop
    L_HX_max : float = 0.545 + 50e-3           # [m] max HX axial length
    mass_HX_max : float = 100          # [kg] max HX mass (optimizer constraint)
    T_c_safety_factor : float = 2.5    # lower T_c limit = T_c_min_coolprop * this factor
    radiation_ON : bool = True         # enable WSGGM thermal radiation (shellnHelicalTube only)
    # Chemistry model: "finite_rate" (FPV default), "equilibrium", or "frozen".
    chemistry_model : str = "finite_rate"
    equilibrium_dh_gas_ON : bool = True  # legacy alias ignored when chemistry_model is explicit
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
class transientProp:
    """Boundary-condition schedules and integration controls for the dynamic solver.

    Each schedule is a list of (t_seconds, value) points, linearly interpolated and
    held flat outside the range; pass None to hold the steady input value constant.
    See DESIGN_PLAN_shellntube_transient.md section 4.  The wall state is a lumped
    thickness-mean temperature per axial node (validated <2 K vs a resolved PDE).
    """
    t_end : float = 40.0                 # [s] total simulated time
    max_step : float = 0.25              # [s] solve_ivp cap; audit tightens it per ramp
    solver_method : str = "BDF"          # "BDF" (implicit, big steps once transient decays) | "RK45"
    n_save : int = 120                   # number of stored time snapshots for the dashboard
    n_axial : int = 80                   # transient axial nodes (coarsened from the steady grid;
                                         #   the wall field is smooth so ~80 resolves it well and
                                         #   keeps the per-RHS march cheap)
    T_wall_initial : float = 293.15      # [K] uniform cold-start wall temperature

    # --- BC schedules: list[(t, value)] or None to hold the steady input ---
    # default scenario = start-up: He ramps 0->full in 2.5 s, gas ignites at t=0 and
    # its convective strength ramps in over 1 s.
    schedule_mass_flow_c = ((0.0, 1e-3), (2.5, 0.15))          # [kg/s] He
    schedule_mass_flow_g = ((0.0, 5e-3), (1.0, 0.10))          # [kg/s] hot gas
    schedule_T_c_in = None               # [K]  hold coolantProp.T_in
    schedule_p_c_in = None               # [Pa] hold coolantProp.p_in
    schedule_OF = None                   # [-]  hold hotgasProp.mixing_ratio
    schedule_T_lox_in = None             # [K] stored for input traceability; not yet used by gas manifold
    schedule_ignition_state = None       # [0/1] stored for input traceability; ignition_time gates physics

    # --- chemistry during transient (see section 4b) ---
    # "equilibrium" (default, REQUIRED for this high-heat-extraction regime — tabulated
    # via the equilibrium manifold so it is as fast as frozen) | "frozen" (validation only)
    # Default transient combustion model: FPV manifold + progress-variable transport.
    chemistry_transient : str = "finite_rate"
    ignition_time : float = 0.0           # [s] gas inert (no flame heat) before this

    # --- cold-side early-ramp caveat flag (section 4.6) ---
    flag_He_outlet_when_residence_gt_tau : bool = True


@dataclass
class runProp:
    """User-facing run controls.

    The intended entry points are:
      python -m hps_combustor.main_steady
      python -m hps_combustor.main_transient

    Select the HX type with combustorProp.HX_config:
      "shellnHelicalTube" or "shellntube".
    """
    run_name: str = "combustor_hx"
    output_root: str = "zip_folders"
    make_archive: bool = True
    save_csv: bool = True
    save_input_snapshot: bool = True
    plot_dashboard: bool = False

    # Optional transient schedule file. CSV is dependency-free; XLSX uses pandas
    # if available. See schedule_inputs.py for accepted columns/sheets.
    schedule_file: Optional[str] = None

    # Shell-and-tube has an independent axial grid from the helical arc grid.
    shelltube_steady_nodes: int = 200
    shelltube_transient_nodes: int = 80


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
    # ---- Shell-and-tube config (WP1) ----
    # Tube-side gas Nu variable-property exponent: Nu *= (T_wall/T_bulk)^n. n=0 default —
    # the hot-end T_b/T_w~2.5 is outside correlation data, so this is a calibration knob.
    n_tube_gas : float = 0.0
    Re_transition_lo : float = 2300.0   # tube-side laminar->transition threshold
    Re_transition_hi : float = 4000.0   # tube-side transition->turbulent threshold
    tube_grooved_Nu_factor : float = 1.0  # EchTherm grooved-tube placeholder; calibrate if using grooves
    tube_grooved_f_factor : float = 1.0   # EchTherm grooved-tube pressure-drop placeholder
    tube_intensification_factor : float = 1.0
    # Bell-Delaware overall prefactor knobs (calibration against shell-side data)
    zukauskas_C_factor : float = 1.0    # overall ideal-bank h prefactor
    bell_Jl_factor : float = 1.0        # leakage-correction multiplier
    bell_Jb_factor : float = 1.0        # bypass-correction multiplier

