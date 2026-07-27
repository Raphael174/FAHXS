""" 
@ author : George
"""
#%%

#!PHYSICS
from physics.combustion_chemistry.combustion_gas import combustion_gas_solve, choose_fuel
from physics.friction_correlations import getFrictionColebrook1939, getFrictionDeveloping, getFrictionCurvedPipeAli2024
from physics.heat_transfer_correlations import getNusseltDeveloping, compute_Nusselt_Gnielinski
from physics.heat_transfer_correlations import nu_churchill_bernstein, phi_sl_multiplier, nusselt_toroid_Ahmed1997, nusselt_inner_curved_tube_mori1967
from physics.heat_conduction import OneDimensionalSteadyConduction_ShellnHelicalTube
from physics.governing_equations import *
from physics.radiation_model.radiation_build import make_ehlme_backend
#!MECHANICAL
from mechanical.material_specs.material_temperature_strength import init_material_temperature_properties

from mechanical.geometry.helix_geometry import *
from mechanical.loads import *
#!DATA
from model_data_process.data_processing import solver_data
from model_data_process.data_plotting import plotALL
#!LIBRARIES
from CoolProp.CoolProp import PropsSI
from CoolProp import AbstractState
import numpy as np
from scipy.constants import gas_constant
from scipy.integrate import simpson
from pathlib import Path

#%%

#%%
class main_solver : 

    def __init__ (  self, 
                    coolantProp,
                    hotgasProp,
                    combustorProp,
                    numericalProp,
                    system_requirements):

        # extract dataclasses 
        self.coolantProp = coolantProp  
        self.hotgasProp = hotgasProp
        self.combustorProp = combustorProp 
        self.numericalProp = numericalProp
        self.system_requirements = system_requirements

        # fuel set up 
        self.chem_mech_path, self.Y_fuel, self.Hv_fuel = choose_fuel(self.hotgasProp.fuel)

        """ 
        DEFINE HELICAL COIL GEOMETRY
        """
        self.N_ch = 1
        self.Dh_ch = self.combustorProp.Dh_coil 
        self.A_ch = np.pi*self.combustorProp.Dh_coil **2/4
        self.coil_pitch = self.Dh_ch + 2*self.combustorProp.thickness_coil_wall + self.combustorProp.coil_gap
        self.D_coil = self.combustorProp.inner_diameter - 2*self.combustorProp.gap_shell2coil - self.combustorProp.Dh_coil - 2*self.combustorProp.thickness_coil_wall
        self.D_inner_coil_passage = self.D_coil - self.combustorProp.Dh_coil - 2*self.combustorProp.thickness_coil_wall
        self.D_tube = self.Dh_ch + 2*self.combustorProp.thickness_coil_wall

        self.Rc = self.D_coil/2*(1 + (self.coil_pitch/(np.pi*self.D_coil))**2)  # coil curvature
        if self.D_coil < 0 : 
            raise Exception("Coil center-to-center diameter is negative - check your geometry")

        # coil axial length from pipe centers
        self.L_coil = (self.numericalProp.L_HX_max-self.combustorProp.mixing_length) - 2*self.combustorProp.length_2_coil - (self.Dh_ch+2*self.combustorProp.thickness_coil_wall)
        #* get functions to translate tube length to axial combustor length and get total tube length
        self.func_s_to_x, self.func_s_to_theta, self.L_ch_max = HelixGeometryRadiusCST(   coil_pitch=self.coil_pitch,
                                                                                D_coil=self.D_coil, 
                                                                                L_coil=self.L_coil)
        self.Dh_cc = compute_Dh_shell(D_coil=self.D_coil, 
                                d_coil_outer=self.Dh_ch+2*self.combustorProp.thickness_coil_wall, 
                                shell_diameter=self.combustorProp.inner_diameter, 
                                coil_pitch=self.coil_pitch)
        self.Ap_cc = np.pi*self.Dh_cc**2/4

        self.area_g_square_PT = self.combustorProp.inner_diameter**2 - np.pi*( (self.combustorProp.inner_diameter - 2*self.combustorProp.gap_shell2coil)**2 - (self.D_coil-self.combustorProp.Dh_coil - 2*self.combustorProp.thickness_coil_wall)**2)/4
        self.area_g_round_PT = np.pi*self.combustorProp.inner_diameter**2/4 - np.pi*( (self.combustorProp.inner_diameter - 2*self.combustorProp.gap_shell2coil)**2 - (self.D_coil-self.combustorProp.Dh_coil - 2*self.combustorProp.thickness_coil_wall)**2)/4
        self.relative_area_dif = (self.area_g_square_PT-self.area_g_round_PT)/self.area_g_round_PT
        print(f'relative diff square vs round PT={self.relative_area_dif}')

        #* HX pipe wall port area for mass computation
        self.Ap_HX = np.pi*((self.Dh_ch+2*self.combustorProp.thickness_coil_wall)**2 - self.Dh_ch**2)/4

        self.P_wc_unit = np.pi*self.combustorProp.Dh_coil 
        #! coil_mass_fraction_g defines mass_fraction of total shell gas on which coil heat transfer acts on
        self.numericalProp.dx = np.pi*self.D_coil
        
        self.coil_mass_fraction_g = self.numericalProp.dx/(np.pi*self.D_coil)

        #! remove restriction on Helium inlet properties
        #* ensures respected geometry of HX

        """
        MODEL ---------------------------------- RADIATION BUILD
        """
        if self.numericalProp.radiation_ON == True : 
            #! RADIATIVE PARAMETERS FOR CO2-H2O MIXTURE 
            self.radiation_backend = make_ehlme_backend(Path(__file__).parent/"physics/radiation_model/ehlme2025_mixture.json")  # if only mixture is present
            #! mean beam length radiative heat
            # taking the approximate effective gas volume around each coil
            V_tot_1_turn = np.pi*self.combustorProp.inner_diameter**2/4 * self.numericalProp.L_HX_max 
            V_pipe = np.pi*(self.Dh_ch+2*self.combustorProp.thickness_coil_wall)**2/4*self.L_ch_max
            A_wg = np.pi*(self.Dh_ch+2*self.combustorProp.thickness_coil_wall)*self.L_ch_max
            self.Le = 3.4*(V_tot_1_turn-V_pipe)/A_wg 
            #!----------------------------------
            print("Le=", self.Le, "V_tot_1_turn-V_pipe=", V_tot_1_turn-V_pipe)
            print(f"Pipe max length = {self.L_ch_max} m")

        else:
            self.radiation_backend = None
        
        
        #! using the shell hydraulic diameter to compute the gas bulk velocity

        #* material properties divided between conductive zone (HX) and burner (CC)
        self.func_CTE_HX, self.func_E_HX, self.func_Yield_HX, self.func_conductivity_HX, self.density_HX, self.poisson_HX, self.func_cp_HX = init_material_temperature_properties(self.combustorProp.material_HX)
        _, _, _, _, self.density_CC, _, _ = init_material_temperature_properties(self.combustorProp.material_CC)


        # initialize coolant #! counter-flow configuration so starting from exit temperature
        self.T_c, self.p_c = self.coolantProp.T_out, self.coolantProp.p_out
        self.rho_c = PropsSI('D','T', self.T_c,'P',self.p_c ,self.coolantProp.coolant)
        self.U_c = self.coolantProp.mass_flow_c / (self.rho_c * self.A_ch * self.N_ch)

        # initialize hot gas 
        self.p_g = np.copy(self.hotgasProp.p0)
        self.combustion_node = combustion_gas_solve(fuel=self.hotgasProp.fuel, oxidizer=self.hotgasProp.oxidizer, OF=self.hotgasProp.mixing_ratio,
                                T_inj_LOX=self.hotgasProp.T_inj_LOX, T_g_init=self.hotgasProp.T_g_init, 
                                p0=self.p_g, 
                                chem_mech_path=self.chem_mech_path, Hv_fuel=self.Hv_fuel, Y_fuel=self.Y_fuel)
        self.combustion_node.solve() # combustion gas properties extractable using cantera
        self.gas_phase = self.combustion_node.phase
        self.T_g  = self.gas_phase.T
        self.rho_g = self.gas_phase.density
        self.U_g = self.hotgasProp.mass_flow_g/(self.rho_g*self.Ap_cc)

        # initialize numerical parameters 
        self.L_HX = 0
        if self.combustorProp.HX_config=="shellnHelicalTube" :
            self.L_HX+=self.Dh_ch+2*self.combustorProp.thickness_coil_wall
        self.L_ch = 0
        self.T_wg, self.T_wc, self.T_c_check = np.copy(self.T_g), np.copy(self.T_c), np.copy(self.T_c) # initialize temperatures for 1D conduction

        # initialize data record
        self.data_master = solver_data

        #! species index for radiation, molar fractions
        self.index_H2O = self.gas_phase.species_index("H2O")
        self.index_CO2 = self.gas_phase.species_index("CO2")

        AS = AbstractState('HEOS',self.coolantProp.coolant)
        self.T_c_min_coolprop = AS.Tmin() # K
        print(f"T_c_min_coolprop= {self.T_c_min_coolprop}")
        print(f"corrected={1.15*self.T_c_min_coolprop}")


    def solver (self) :

        # iterate through HX using numerical parameters 
        #* stay in loop as long as :
            #! HX below max length
            #! T_c above CoolProp min limit
        while   self.L_HX  <= (self.numericalProp.L_HX_max - self.combustorProp.mixing_length - 2*self.combustorProp.length_2_coil) \
                and self.T_c > self.T_c_min_coolprop*2.5:

            """ 
            Resolve heat transfer 
            """

            #* check for compressibility factor
            self.Z = PropsSI('Z','T',self.T_c,'P',self.p_c,self.coolantProp.coolant)
            # if abs(self.Z-1) > 0.03:  # 3%
            #     print(f"Warning - Z={self.Z}")            # start with initialized (U, p, T, rho)_c and T_g
            
            # Coolant thermodynamics + convection coef
            self.cp_c = PropsSI('C','T', self.T_c,'P',self.p_c,self.coolantProp.coolant)
            self.cv_c = PropsSI('CVMASS','T', self.T_c,'P',self.p_c,self.coolantProp.coolant)
            self.gamma_c = self.cp_c/self.cv_c
            self.mu_c = PropsSI('V','T',self.T_c,'P',self.p_c,self.coolantProp.coolant)
            self.k_c = PropsSI('L','T',self.T_c,'P',self.p_c,self.coolantProp.coolant)
            # flow characteristics
            self.Re_c=  self.rho_c*self.U_c*self.Dh_ch/self.mu_c
            self.De = self.Re_c*np.sqrt(self.Dh_ch/self.D_coil) # Dean number
            self.He = self.Re_c*np.sqrt(self.Dh_ch/(2*self.Rc)) # Helical number
            """ 
            Coil tube friction correlation
            """
            if self.combustorProp.friction_coil == "Colebrook1939":
                self.f_fd_c = getFrictionColebrook1939(self.Re_c, self.combustorProp.channel_roughness/self.Dh_ch)
            elif self.combustorProp.friction_coil == "CurvedPipeAli2024":
                self.f_fd_c = getFrictionCurvedPipeAli2024(Re=self.Re_c, Dh=self.Dh_ch, Rc=self.D_coil/2)
            else:
                self.f_fd_c = getFrictionColebrook1939(self.Re_c, self.combustorProp.channel_roughness/self.Dh_ch)
                raise ('Warning - coil friction correlation not found - using Colebrook')
            # check for fully rough region 
            self.f_c = getFrictionDeveloping(self.f_fd_c, self.Dh_ch, 10e10)*(1 + self.numericalProp.artificial_error_friction_cold)
            
            """ 
            Coil Nusselt correlations
            """
            # coolant heat transfer properties
            self.Pr_c = self.cp_c*self.mu_c/self.k_c
            if self.combustorProp.Nusselt_coil == "Gnielinski":
                self.Nu_c = compute_Nusselt_Gnielinski(f= self.f_fd_c, Re_c= self.Re_c, Pr_c=self.Pr_c)
            elif self.combustorProp.Nusselt_coil == "mori1967":
                self.Nu_c = nusselt_inner_curved_tube_mori1967(Re=self.Re_c, Pr=self.Pr_c, d=self.Dh_ch, R=self.D_coil/2)
            else:
                self.Nu_c = compute_Nusselt_Gnielinski(f= self.f_fd_c, Re_c= self.Re_c, Pr_c=self.Pr_c)
                raise ('Warning - coil Nusselt correlation not found - using Gnielinski')
            # compute HT
            self.h_c = getNusseltDeveloping(self.Nu_c, self.Dh_ch, 10e10)*self.k_c/self.Dh_ch

            """ 
            Hot gas Nusselt correlations
            """
            self.cp_g, self.cv_g = self.gas_phase.cp, self.gas_phase.cv
            self.gamma_g = self.cp_g/self.cv_g
            self.mu_g = self.gas_phase.viscosity 
            self.k_g = self.gas_phase.thermal_conductivity
            self.rho_g = self.gas_phase.density
            self.W_g = self.gas_phase.mean_molecular_weight
            # flow characteristics
            self.U_g = self.hotgasProp.mass_flow_g/(self.rho_g*self.Ap_cc)
            self.Re_g= self.rho_g*self.U_g*self.Dh_cc/self.mu_g
            self.f_fd_g = getFrictionColebrook1939(self.Re_g, self.combustorProp.combustor_roughness/self.Dh_cc)
            # check for fully rough region 
            self.f_g = getFrictionDeveloping(self.f_fd_g, self.Dh_cc, 10e10)
            # coolant heat transfer properties
            self.Pr_g = self.cp_g*self.mu_g/self.k_g

            #! radiation inputs
            self.X_H2O = self.combustion_node.phase.X[self.index_H2O]
            self.X_CO2 = self.combustion_node.phase.X[self.index_CO2]


            self.Re_sh = self.rho_g*self.U_g*(self.Dh_ch+2*self.combustorProp.thickness_coil_wall)/self.mu_g

            if self.combustorProp.Nusselt_shell == "churchill_bernstein":
                self.Nu_g = nu_churchill_bernstein(Re=self.Re_sh, Pr=self.Pr_g)
                self.h_g = self.Nu_g*self.k_g/(self.Dh_ch+2*self.combustorProp.thickness_coil_wall)

            elif self.combustorProp.Nusselt_shell == "churchill_bernstein_tightcoil":
                self.phi_multiplier = phi_sl_multiplier(SL_over_D=self.coil_pitch/(self.Dh_ch+2*self.combustorProp.thickness_coil_wall))
                self.Nusselt_correction = self.combustorProp.Nusselt_correction
                self.Nu_g = self.Nu_g*self.phi_multiplier*self.Nusselt_correction
                
                #! RE WRITE THIS WITH PROPER FORM -THEN TIGHT COIL BECOMES UPPER BOUND
                
                self.h_g = self.Nu_g*self.k_g/(self.Dh_ch+2*self.combustorProp.thickness_coil_wall)
            
            elif    self.combustorProp.Nusselt_shell == "ahmed_toroid" or \
                    self.combustorProp.Nusselt_shell != "churchill_bernstein" and \
                    self.combustorProp.Nusselt_shell != "churchill_bernstein_tightcoil"  :
                self.Asqrt_toroid = np.sqrt((self.Dh_ch+2*self.combustorProp.thickness_coil_wall)*self.D_coil*np.pi**2) #perimeter tube x perimeter coil
                self.Nu_g = nusselt_toroid_Ahmed1997 (U_g=self.U_g, rho_g=self.rho_g, mu_g=self.mu_g, Pr_g=self.Pr_g, Asqrt_toroid=self.Asqrt_toroid)
                self.h_g = self.Nu_g*self.k_g/self.Asqrt_toroid
            


            self.c_c = np.sqrt(gas_constant*1e3/self.coolantProp.molar_mass*self.T_c*self.cp_c/self.cv_c)
            self.Mach_c = self.U_c/self.c_c

            self.c_g = np.sqrt(gas_constant*1e3/self.W_g*self.T_g*self.cp_g/self.cv_g)
            self.Mach_g = self.U_g/self.c_g
            # 1D heat conduction 

            self.heat_transfer_node =  OneDimensionalSteadyConduction_ShellnHelicalTube(
                                        h_g=self.h_g*(1 + self.numericalProp.artificial_error_Nu_hot), h_c=self.h_c*(1 + self.numericalProp.artificial_error_Nu_cold),
                                        T_c=self.T_c, T_g=self.T_g, 
                                        s_w=self.combustorProp.thickness_coil_wall,
                                        Dh_ch=self.Dh_ch,
                                        f_kw_at_T=self.func_conductivity_HX,
                                        T_wg_0=self.T_wg, T_wc_0=self.T_wc, T_c_check_0=self.T_c_check,
                                        dx=self.numericalProp.dx,
                                    # radiation:
                                        rad_enabled=self.numericalProp.radiation_ON,
                                        eps_s=self.numericalProp.emissivity_wall,
                                        rad_backend=self.radiation_backend,
                                        rad_state={'p': self.p_g, 'yH2O': self.X_H2O, 'yCO2': self.X_CO2, 'Le': self.Le})
            

            self.heat_transfer_node.Solve1Dconduction()

            self.UP = self.heat_transfer_node.UP
            self.UA = self.heat_transfer_node.UA
            self.dQ = self.heat_transfer_node.dQ #! heat extracte from this cell, used to extract heat from hot gases
            self.q_w = self.heat_transfer_node.q_w #self.dQ/(self.numericalProp.dx*self.N_ch*(self.w_ch+2*self.h_ch)) # heat flux according to nodal power on channel conducting surfaces
            self.dq__dx = self.heat_transfer_node.dq__dx #! +ive here, but counter flow means reversed gradient direction
            self.Res_g = self.heat_transfer_node.Res_g 
            self.Res_c = self.heat_transfer_node.Res_c 
            self.Res_w = self.heat_transfer_node.Res_w 
            self.Biot_g = self.Res_w/self.Res_g
            self.Biot_c = self.Res_w/self.Res_c
            #! radiation data [CODED FOR HELICAL TUBE HX]
            self.q_w_rad = self.heat_transfer_node.q_w_rad
            self.emissivity_g = self.heat_transfer_node.eps_emit 
            self.absorptivity_g = self.heat_transfer_node.eps_abs
            self.h_g_rad = self.heat_transfer_node.h_g_rad
            self.h_g_conv = self.heat_transfer_node.h_g
            # print("q_w_rad=", self.q_w_rad,"h_g_rad=", self.h_g_rad , "h_g_conv=", self.h_g_conv)
            # print(f"emissivity_g={self.emissivity_g}, absorptivity_g={self.absorptivity_g}")
            # print(f"T_g={self.T_g}, T_wg={self.T_wg}")
            if self.combustorProp.HX_config == "shellnHelicalTube":
                self.dh_g = self.dQ/(self.hotgasProp.mass_flow_g*self.coil_mass_fraction_g)
            else: 
                self.dh_g = self.dQ/self.hotgasProp.mass_flow_g
            
            # wall temperatures
            self.T_wg = self.heat_transfer_node.T_wg_new
            self.T_wc = self.heat_transfer_node.T_wc_new
            self.k_w = self.heat_transfer_node.k_w
            # check cold fluid temperatute 
            self.T_c_check = self.heat_transfer_node.T_c_check_f # used to double check resolution
            
            """ 
            STRESSES
            """
            #!extract functions
            #* properties at mean wall temperatures
            self.CTE = self.func_CTE_HX((self.T_wg+self.T_wc)/2-273)
            self.Modulus = self.func_E_HX((self.T_wg+self.T_wc)/2-273)
            self.Yield = self.func_Yield_HX((self.T_wg+self.T_wc)/2-273)
            #* stresses
            self.stress_pressure = stress_pressure_tube(P=self.p_c, thickness_pipe=self.combustorProp.thickness_coil_wall, Dh_pipe=self.Dh_ch)
            self.stress_thermal_inner,  self.stress_thermal_outer= stress_thermal_tube(T_inner=self.T_wc, T_outer=self.T_wg, 
                                                                                        CTE=self.CTE,
                                                                                        E=self.Modulus, 
                                                                                        poisson=self.poisson_HX)
            self.stress_inner, self.stress_outer = [self.stress_thermal_inner+self.stress_pressure,  self.stress_thermal_outer+self.stress_pressure]

            # compute derivatives of coolant 
            #* ignore impact of pressure on temperature for now 
            #! using power gradient dq__dx in order to conserve the effective perimeters (i.e. heat exchange widths)
            self.dU_c__dx = dU__dx_IdealGas(U=self.U_c, A=self.A_ch, p=self.p_c, P_w=1, q_w=self.dq__dx/self.N_ch, T=self.T_c, cp=self.cp_c, m_dot=self.coolantProp.mass_flow_c/self.N_ch, dA__dx=0, f=self.f_c, Dh=self.Dh_ch)
            self.dT_c__dx = dT__dx_IdealGas(P_w=1, q_w=self.dq__dx/self.N_ch, m_dot=self.coolantProp.mass_flow_c/self.N_ch, U=self.U_c, dU__dx=self.dU_c__dx, cp=self.cp_c)
            #self.dp_c__dx = dp__dx_IdealGas(p=self.p_c, T=self.T_c, dT__dx=self.dT_c__dx, U=self.U_c, dU__dx=self.dU_c__dx, A=self.A_ch, dA__dx=0)
            #self.drho_c__dx = drho__dx_IdealGas(rho=self.rho_c, p=self.p_c, dp__dx=self.dp_c__dx, T=self.T_c, dT__dx=self.dT_c__dx)
            self.drho_c__dx = drho__dx_IdealGas_logical(rho=self.rho_c, U=self.U_c, dU__dx=self.dU_c__dx, A=self.A_ch, dA__dx=0)
            self.dp_c__dx = dp__dx_IdealGas_logical(p=self.p_c, T=self.T_c, dT__dx=self.dT_c__dx, rho=self.rho_c, drho__dx=self.drho_c__dx)
            # self.dT_c__dx = -(self.T_g - self.T_c)*self.UP/(self.coolantProp.mass_flow_c*self.cp_c)
            # self.dp_c__dx = self.f_c*self.rho_c*self.U_c**2/(2*self.Dh_cc)
            # compute dp hot gas 
            self.dp_g__dx = -self.f_g*self.rho_g*self.U_g**2/(2*self.Dh_cc)

            #print("HX_mass = ", self.HX_mass, " kg")
            #! record data from master dict
            for key in self.data_master:
                #* check if key exists as attribute otherwise skip
                if hasattr(self, key):
                    self.data_master[key].append(getattr(self, key))

            #! update state variables after data record

            """
            Track lengths
            """
            self.L_ch+=self.numericalProp.dx 
            if self.combustorProp.HX_config == "shellnHelicalTube" or self.combustorProp.HX_config == "coolingcoil":
                self.L_HX = self.func_s_to_x(self.L_ch) + self.Dh_ch+2*self.combustorProp.thickness_coil_wall
            else:
                self.L_HX+=self.numericalProp.dx 

            """
            Update state variables coolant
            """
            self.T_c+= self.dT_c__dx*self.numericalProp.dx *-1
            #*self.rho_c+= self.drho_c__dx*self.numericalProp.dx 
            self.p_c+= self.dp_c__dx*self.numericalProp.dx  *-1
            self.rho_c+= self.drho_c__dx*self.numericalProp.dx  *-1
            # self.rho_c = PropsSI('D','T', self.T_c,'P', self.p_c ,self.coolantProp.coolant)
            self.U_c+= self.dU_c__dx*self.numericalProp.dx  *-1
            # self.U_c = self.coolantProp.mass_flow_c/(self.rho_c*self.A_ch*self.N_ch)
            """
            Update state combustion gas 
            """
            self.p_g+=self.dp_g__dx*self.numericalProp.dx 
            self.combustion_node.remove_energy(dh=self.dh_g, updated_pressure=self.p_g, equilibrium_dh_gas_ON=self.numericalProp.equilibrium_dh_gas_ON)
            self.T_g = self.combustion_node.phase.T 
            self.cp_g, self.cv_g = self.combustion_node.phase.cp, self.combustion_node.phase.cv
            self.mu_g = self.combustion_node.phase.viscosity 
            self.k_g = self.combustion_node.phase.thermal_conductivity
            self.rho_g = self.combustion_node.phase.density
            


    
    def HX_sizing_brief (self, plotON=True, printON=True):
        
        if plotON==True :
            plotALL(self.data_master)
                
        """
        HX PERFORMANCE
        """
        self.C_g_avg = np.average(np.array(self.data_master["cp_g"])*self.hotgasProp.mass_flow_g)
        self.C_c_avg = np.average(np.array(self.data_master["cp_c"])*self.coolantProp.mass_flow_c)
        

        self.T_wg_max = np.max(self.data_master["T_wg"])
        self.T_wc_max = np.max(self.data_master["T_wc"])

        self.dT_max = np.max(np.array(self.data_master["T_wg"]) - np.array(self.data_master["T_wc"]))

        self.eta_HX = self.C_g_avg*(self.data_master["T_g"][0]-self.data_master["T_g"][-1])/(np.min([self.C_c_avg, self.C_g_avg])*(self.data_master["T_g"][0]-self.data_master["T_c"][-1]))
        self.Mach_g_max = np.max(self.data_master["Mach_g"])
        self.Q_tot = sum(self.data_master["dQ"])*1e-3
        self.Q_He = np.abs(self.coolantProp.mass_flow_c * simpson(np.array(self.data_master["cp_c"]), np.array(self.data_master["T_c"])))*1e-3
        self.T_fin_max = np.max(self.data_master["T_g"])
        self.dp_c_tot = self.data_master["p_c"][-1]-self.data_master["p_c"][0]
        
        """ 
        check for combustion gases propulsion effects
        """
        self.PR_crit_g = (2/(self.gamma_g+1))**(self.gamma_g/(self.gamma_g-1)) # critical downstream pressure/upstream total pressure
        self.p0_choke = self.system_requirements.ambient_pressure/self.PR_crit_g
        self.Gam_g = np.sqrt(self.gamma_g)*(2/(self.gamma_g+1))**((self.gamma_g+1)/(2*(self.gamma_g-1)))
        self.A_sonic_g = self.hotgasProp.mass_flow_g*np.sqrt(gas_constant*1e3/self.W_g*self.T_g)/(self.Gam_g*self.p0_choke)
        self.D_sonic_g = np.sqrt(4/np.pi*self.A_sonic_g)
        self.D_pipe_g = self.combustorProp.exhaust_diameter
        self.v_exit_g = self.hotgasProp.mass_flow_g/(self.rho_g*np.pi*self.D_pipe_g**2/4)
        self.Thrust_burner = self.v_exit_g*self.hotgasProp.mass_flow_g 
        self.v_exit_g_required = self.system_requirements.max_thrust/self.hotgasProp.mass_flow_g 
        self.A_exit_required = self.hotgasProp.mass_flow_g/(self.v_exit_g_required*self.rho_g)
        self.D_exit_required = np.sqrt(4/np.pi*self.A_exit_required)

    

        # mass propellant 
        self.F = 1/(self.hotgasProp.mixing_ratio + 1)
        self.O = 1 - self.F
        self.mass_kerosene = self.system_requirements.burn_time*self.hotgasProp.mass_flow_g * self.F
        self.mass_LOx = self.system_requirements.burn_time*self.hotgasProp.mass_flow_g * self.O

        #* masses identiques pour toutes solutions
        # injector plate i.e. inlet cap
        self.mass_injector_plate = self.density_CC*(np.pi*self.combustorProp.inner_diameter**2/4)*self.combustorProp.wall_thickness_inj
        # mass gas ejection, x1.1 the surface area of the cylinder cross section
        self.mass_gas_eject = self.density_CC*(1.1*np.pi*self.combustorProp.inner_diameter**2/4)*self.combustorProp.wall_thickness_cc
        # combustor wall mass mixing length
        self.mass_mixing_zone = self.density_CC*self.combustorProp.mixing_length*np.pi*((self.combustorProp.inner_diameter+2*self.combustorProp.wall_thickness_cc)**2 - self.combustorProp.inner_diameter**2)/4
        #* specific to each
        self.N_turns = self.L_coil/self.coil_pitch
        #* HX pipe mass
        self.mass_HX = self.L_ch * self.density_HX * np.pi*((self.Dh_ch+2*self.combustorProp.thickness_coil_wall)**2 - self.Dh_ch**2)/4
        #* HX shell mass only around HX zone
        self.mass_shell_walls = self.density_CC*(self.L_HX-self.combustorProp.mixing_length)*np.pi*((self.combustorProp.inner_diameter + 2*self.combustorProp.wall_thickness_cc)**2 - self.combustorProp.inner_diameter**2)/4
        self.mass_combustor = self.mass_injector_plate + self.mass_gas_eject + self.mass_mixing_zone + self.mass_shell_walls

        self.mass_tot = self.mass_kerosene + self.mass_LOx + self.mass_HX + self.mass_combustor 

        """
        Record Stresses
        """
        # compute max stress/yield ratio, this max should be kep low, say 50%
        self.max_stress__yield = np.max([np.max( np.array(self.data_master["stress_inner"])/np.array(self.data_master["Yield"])), np.max( np.array(self.data_master["stress_outer"])/np.array(self.data_master["Yield"]))] )
        
        if printON==True:
            print(" ")
            print("RESULTS : ")
            print("HX configuration :", self.combustorProp.HX_config)
            print("N_ch = ", self.N_ch)
            if self.combustorProp.HX_config == "coolingjacket" or self.combustorProp.HX_config == "shellntube":
                print("mL_cc = ", self.heat_transfer_node.m_cc*self.combustorProp.fin_height)
            print("L_pipe = ", self.L_ch)
            print("L_HX = ", self.L_HX)
            print("D_coil = ", self.D_coil)
            print("N_turns = ", self.L_ch/(np.pi*self.D_coil))
            print(" ")    
            print("eta_HX =", self.eta_HX )
            print("T_c_in = ", self.T_c, " K")
            print("p_c_in = ", self.p_c/1e5, " bar")
            print("coil curvature = ", self.Rc)
            print(" ")
            print("D_inner_coil_passage = ", self.D_inner_coil_passage*1e3, " mm")
            print("D_tube = ", self.D_tube*1e3, " mm")
            print("D_inner_coil_passage/D_tube = ", (self.D_inner_coil_passage/self.D_tube), " mm")
            print("  ")
            print("D_coil/Dh= ", self.D_coil/self.Dh_ch)
            print("Dean_max = ", np.max(self.data_master["De"]), " Dean_min = ", np.min(self.data_master["De"]))
            print("Helical_number_max = ", np.max(self.data_master["He"]), " Helical_number_min = ", np.min(self.data_master["He"]))
            print("Nu_c_max = ", np.max(self.data_master["Nu_c"]), " Nu_c_min = ", np.min(self.data_master["Nu_c"]))
            print("Nu_g_max = ", np.max(self.data_master["Nu_g"]), " Nu_g_min = ", np.min(self.data_master["Nu_g"]))
            print("Biot_g_max = ", np.max(self.data_master["Biot_g"]), " Biot_g_min = ", np.min(self.data_master["Biot_g"]))
            print("Biot_c_max = ", np.max(self.data_master["Biot_c"]), " Biot_c_min = ", np.min(self.data_master["Biot_c"]))
            if self.combustorProp.Nusselt_shell == "churchill_bernstein_tightcoil": 
                print("phi_max = ", np.max(self.data_master["phi_multiplier"]), " phi_min = ", np.min(self.data_master["phi_multiplier"]))
            print("Re_c_max = ", np.max(self.data_master["Re_c"])/1e6, "e6")
            print("dP = ", (self.p_c-self.coolantProp.p_out)/1e5, " bar")
            print("L_HX = ", self.L_HX, " m")
            print("Q_tot = ", self.Q_tot, "Q_He = ", self.coolantProp.mass_flow_c*np.average(self.data_master["cp_c"])*(self.data_master["T_c"][0]-self.data_master["T_c"][-1])/1e3, self.Q_He, " kw")
            print("T_wg_max = ", self.T_wg_max, " K")
            print("T_wc_max = ", self.T_wc_max, " K")
            print("dT_wall_max = ", self.dT_max)
            print("T_g_max = ", np.max(self.data_master["T_g"]), " K")
            print("Mach_g max = ", self.Mach_g_max)
            print(" ")
            print("mass combustor-HX = ", self.mass_combustor+self.mass_HX, " kg")
            print("mass Liquid Fuel = ", self.mass_kerosene, " kg")
            print("mass LOX = ", self.mass_LOx, " kg")
            print("total = ", self.mass_tot, " kg")
            print(" ")
            print("Combustion gases choking conditions")
            print(f"Choking pressure = {self.p0_choke/1e5} bar @ p_amb={self.system_requirements.ambient_pressure/1e5} bar")
            print(f"Critical choking hot gas diameter = {self.D_sonic_g*1e3} mm")
            print(f"V_exit_g = {self.v_exit_g} m/s | thrust = {self.Thrust_burner} N // assuming p0=ambient")
            print(f"To maintain {self.system_requirements.max_thrust}N, Exit diameter >= {self.D_exit_required*1e3}mm")
            
            print("  ")
            print(" STRESSES ")
            print("max_inner_stress = ", np.max(self.data_master["stress_inner"])/1e6, " MPa")
            print("max_outer_stress = ", np.max(self.data_master["stress_outer"])/1e6, " MPa")
            print(" ")
    
#%%
if __name__ == "__main__":

    from input_data import *

    combustor = main_solver(hotgasProp=hotgasProp, 
                            coolantProp=coolantProp,
                            combustorProp=combustorProp,
                            system_requirements=system_requirements,
                            numericalProp=numericalProp)
    
    combustor.solver()
    combustor.HX_sizing_brief(False)
# %%
