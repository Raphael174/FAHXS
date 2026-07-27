""" 
@ author : Raphaël Aubry
"""

""" 
Script to estimate combustion chemistry using cantera
"""
#%%
import cantera as ct 
from pathlib import Path
import numpy as np 

#from input_data import hotgasProp
#ct.suppress_thermo_warnings()

#%%

base_path = Path(__file__).parent


def choose_fuel(fuel):

    if fuel == "POSF10325" :
        chem_mech_path = base_path / "A2highT.yaml"
        Hv_fuel = 360e3 #J/kg  https://web.stanford.edu/group/haiwanglab/HyChem/approach/Report_Jet_Fuel_Thermochemical_Properties_v6.pdf
        Y_fuel = {'POSF10325' : 1.}

    elif fuel == "H2" :
        chem_mech_path = base_path / "H2-O2_Burke2012.yaml"
        Hv_fuel = 0
        Y_fuel = {'H2' : 1.}

    elif fuel == "gasoline-E5" :
        chem_mech_path = base_path / "llnl_gasoline_323.yaml"
        Y_fuel = {'C2H5OH': np.float64(0.022521123344674744), 'C6H12-1': np.float64(0.05760037978051016),
                'C5H10-1': np.float64(0.020571564207325047),'IC4H8': np.float64(0.010971500910573358),
                'IC8H18': np.float64(0.6121661440853741), 'NC7H16': np.float64(0.2761692876715427)}
        Hv_fuel = 321.219e3 # J/kg

    elif fuel == "gasoline-E10" :
        chem_mech_path = base_path / "llnl_gasoline_323.yaml"
        Y_fuel = {'C2H5OH': np.float64(0.04648289377205732), 'C6H12-1': np.float64(0.059442690615179886),
                'C5H10-1': np.float64(0.021229532362564234), 'IC4H8': np.float64(0.01132241726003426),
                'IC8H18': np.float64(0.5936889008300161), 'NC7H16': np.float64(0.2678335651601482)}
        Hv_fuel = 313.664e3 # J/kg

    elif fuel == "diesel-C16H34" :
        # based on 10.1016/j.fuel.2017.07.009
        # DC surrogate, 
        chem_mech_path = base_path / "RenKokjohn_surrogate.yaml"
        Y_fuel = {'c16h34': np.float64(1)}
        Hv_fuel = 350e3 # J/kg, 10.4271/2008-01-1379
        

    else:
        raise("Error : fuel species not supported")
    
    return chem_mech_path, Y_fuel, Hv_fuel

#%%

class combustion_gas_solve :

    def __init__ (self, 
                  fuel, oxidizer,
                  OF,
                  p0,
                  T_g_init,
                  T_inj_LOX,
                  chem_mech_path, Y_fuel, Hv_fuel):
        
        self.fuel = fuel 
        self.oxidizer = oxidizer
        self.OF = OF 
        self.p0 = p0 
        self.T_g_init = T_g_init 
        self.T_inj_LOX = T_inj_LOX 

        self.phase = ct.Solution(chem_mech_path)
        self.Hv_fuel = Hv_fuel
        self.Y_fuel = Y_fuel

    def solve (self) :

        O = self.OF/(1 + self.OF)
        F = 1 - O

        """ENERGY LOSS LOX-->GOX"""

        # LOX TO GOX
        O2 = ct.Oxygen() # import oxygen 
        O2.TP = self.T_inj_LOX, self.p0 # set O2 to injection thermo
        self.h_ox_inj = O2.enthalpy_mass # fetch O2 onjection enthalpy
        O2.TP = self.T_g_init, self.p0 # set O2 to reaction state, assumed at T_pyro
        self.h_ox_pyro = O2.enthalpy_mass
        self.dH_ox = self.h_ox_pyro-self.h_ox_inj

        self.dH_tot = O*self.dH_ox + F*self.Hv_fuel


        """CHEMICAL EQUILIBRIUM CALCULATION"""

        if self.fuel == "gasoline-E5" or self.fuel == "gasoline-E10" or self.fuel == "diesel-C16H34":
            #! building the O2/gasoline mixture
            s = sum(self.Y_fuel.values())
            if s <= 0: raise ValueError("Y_fuel is empty.")
            Yf = {k:v/s for k,v in self.Y_fuel.items()}
            Y_mix = {k: (1.0/(1.0+self.OF))*Yf[k] for k in Yf}
            Y_mix[self.oxidizer] = self.OF/(1.0+self.OF)

        else:
            #! building the pure fuel/O2 mixture
            Y_mix = "O2" + ':' + str(self.OF) + ',' +  self.fuel + ': 1' 

        self.phase.TPY = self.T_g_init, self.p0, Y_mix

        self.phase.equilibrate('HP') # finding chem equilibrium

        self.phase.HP = self.phase.enthalpy_mass-self.dH_tot, self.p0 # remove energy loss
        self.phase.equilibrate('HP') # re-equilibriate mixture

    equilibrium_dh_gas_ON : bool = True
    def remove_energy (self, dh, updated_pressure, equilibrium_dh_gas_ON): 

        # enthalpy to remove
        self.dh = dh 
        self.p0=updated_pressure

        # equilibriate with removed energy
        self.phase.HP = self.phase.enthalpy_mass-self.dh, self.p0 # remove energy loss
        if equilibrium_dh_gas_ON==True:
            self.phase.equilibrate('HP') # re-equilibriate mixture

#%%

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import pandas as pd

    OF_study = False 
    dQ_study = True

    if OF_study == True :

        fuel = "gasoline-E10"


        chem_mech_path, Y_fuel, Hv_fuel = choose_fuel(fuel)

            # oxygen injection pressure and temperature
        T_inj = 300 # K

        pcomb = 5e5

        T_ig = 800+273.15


        OF_list = np.linspace(2, 4, 50, endpoint=True)
        T_g = np.zeros(len(OF_list))
        cp_g = np.zeros(len(OF_list))
        lambda_g = np.zeros(len(OF_list))
        density_g = np.zeros(len(OF_list))

        for i in range(len(OF_list)):

            
            combustion_object = combustion_gas_solve(
                        fuel=fuel, oxidizer="O2",
                        OF=OF_list[i],
                        p0=pcomb,
                        T_g_init=T_ig,
                        T_inj_LOX=T_inj,
                        chem_mech_path=chem_mech_path, Hv_fuel=Hv_fuel, Y_fuel=Y_fuel)
            combustion_object.solve()
            combustion_object.phase.T
            T_g[i] = combustion_object.phase.T

        OF_max_Tg = OF_list[np.argmax(T_g)]
        # OF_st = 3.42

        plt.plot(OF_list, T_g)
        plt.ylabel(r"$T_g$ [K]")
        plt.xlabel("O/F")
        plt.xlim(OF_list[0], OF_list[-1])
        plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
        # Show the minor grid as well. Style it in very light gray as a thin,
        plt.axvline(x=OF_max_Tg, color='orange', linestyle='--', linewidth=1, label=f"OF@T_max={np.round(OF_max_Tg, 3)}")
        # plt.axvline(x=OF_st, color='blue', linestyle='--', linewidth=1, label="OF@stoich")
        # dotted line.
        plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
        # Make the minor ticks and gridlines show.
        plt.minorticks_on()
        plt.legend()
        plt.show()

    if dQ_study == True :
        """ 
        Study to compute equilibrium and non-equilibrium gas proprties of hot gases after removing removed_power/mass_flow enthalpy 
        """
        # power removed constant
        removed_power = 441504
        # mass flow range for study
        mass_flow_range = np.linspace(75, 150, 50, endpoint=True)*1e-3
        # resulting enthalpy
        dH_range = removed_power/mass_flow_range 
        # thermodynamic metrics
        T_g_equil, T_g_neq = np.zeros(len(dH_range)), np.zeros(len(dH_range))
        cp_g_equil, cp_g_neq = np.zeros(len(dH_range)), np.zeros(len(dH_range))
        cv_g_equil, cv_g_neq = np.zeros(len(dH_range)), np.zeros(len(dH_range))
        lambda_g_equil, lambda_g_neq = np.zeros(len(dH_range)), np.zeros(len(dH_range))
        density_g_equil, density_g_neq = np.zeros(len(dH_range)), np.zeros(len(dH_range))
        viscosity_g_equil, viscosity_g_neq = np.zeros(len(dH_range)), np.zeros(len(dH_range))
        MW_g_equil, MW_g_neq = np.zeros(len(dH_range)), np.zeros(len(dH_range))
        
        # combustion set up 
        fuel = "diesel-C16H34"
        chem_mech_path, Y_fuel, Hv_fuel = choose_fuel(fuel)
        # oxygen injection pressure and temperature
        T_inj = 300 # K
        pcomb = 5e5
        T_ig = 800+273.15
        O__F = 2.9

        equilibrium_dh_gas = True

        for i in range(len(dH_range)) :

            combustion_object = combustion_gas_solve(
                        fuel=fuel, oxidizer="O2",
                        OF=O__F,
                        p0=pcomb,
                        T_g_init=T_ig,
                        T_inj_LOX=T_inj,
                        chem_mech_path=chem_mech_path, Hv_fuel=Hv_fuel, Y_fuel=Y_fuel)
            combustion_object.solve()
            combustion_object.remove_energy(dh=dH_range[i], updated_pressure=pcomb, equilibrium_dh_gas_ON=equilibrium_dh_gas)

            T_g_equil[i] = combustion_object.phase.T
            cp_g_equil[i] = combustion_object.phase.cp 
            cv_g_equil[i] = combustion_object.phase.cv
            lambda_g_equil[i] = combustion_object.phase.thermal_conductivity
            density_g_equil[i] = combustion_object.phase.density
            viscosity_g_equil[i] = combustion_object.phase.viscosity
            MW_g_equil[i] = combustion_object.phase.mean_molecular_weight

        equilibrium_dh_gas = False

        for j in range(len(dH_range)) :

            combustion_object = combustion_gas_solve(
                        fuel=fuel, oxidizer="O2",
                        OF=O__F,
                        p0=pcomb,
                        T_g_init=T_ig,
                        T_inj_LOX=T_inj,
                        chem_mech_path=chem_mech_path, Hv_fuel=Hv_fuel, Y_fuel=Y_fuel)
            combustion_object.solve()
            combustion_object.remove_energy(dh=dH_range[j], updated_pressure=pcomb, equilibrium_dh_gas_ON=equilibrium_dh_gas)


            T_g_neq[j] = combustion_object.phase.T
            cp_g_neq[j] = combustion_object.phase.cp 
            cv_g_neq[j] = combustion_object.phase.cv
            lambda_g_neq[j] = combustion_object.phase.thermal_conductivity
            density_g_neq[j] = combustion_object.phase.density
            viscosity_g_neq[j] = combustion_object.phase.viscosity
            MW_g_neq[j] = combustion_object.phase.mean_molecular_weight
#%%

        #"""TEMPERATURE"""
        plt.plot(mass_flow_range*1e3, T_g_equil, label = "T_g_equil")
        plt.plot(mass_flow_range*1e3, T_g_neq, label = "T_g_neq")
        plt.ylabel(r"$T_g$ [K]")
        plt.xlabel("mass flow [g/s]")
        plt.xlim(mass_flow_range[0]*1e3, mass_flow_range[-1]*1e3)
        plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
        plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
        # Make the minor ticks and gridlines show.
        plt.minorticks_on()
        plt.legend()
        plt.show()
        #"""CP"""
        plt.plot(mass_flow_range*1e3, cp_g_equil, label = "cp_g_equil")
        plt.plot(mass_flow_range*1e3, cp_g_neq, label = "cp_g_neq")
        plt.ylabel(r"$cp_g$ [K]")
        plt.xlabel("mass flow [g/s]")
        plt.xlim(mass_flow_range[0]*1e3, mass_flow_range[-1]*1e3)
        plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
        plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
        # Make the minor ticks and gridlines show.
        plt.minorticks_on()
        plt.legend()
        plt.show()
        #"""CV"""
        plt.plot(mass_flow_range*1e3, cv_g_equil, label = "cv_g_equil")
        plt.plot(mass_flow_range*1e3, cv_g_neq, label = "cv_g_neq")
        plt.ylabel(r"$cv_g$ [K]")
        plt.xlabel("mass flow [g/s]")
        plt.xlim(mass_flow_range[0]*1e3, mass_flow_range[-1]*1e3)
        plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
        plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
        # Make the minor ticks and gridlines show.
        plt.minorticks_on()
        plt.legend()
        plt.show()
        #"""LAMBDA"""
        plt.plot(mass_flow_range*1e3, lambda_g_equil, label = "lambda_g_equil")
        plt.plot(mass_flow_range*1e3, lambda_g_neq, label = "lambda_g_neq")
        plt.ylabel(r"$lambda_g$ [K]")
        plt.xlabel("mass flow [g/s]")
        plt.xlim(mass_flow_range[0]*1e3, mass_flow_range[-1]*1e3)
        plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
        plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
        # Make the minor ticks and gridlines show.
        plt.minorticks_on()
        plt.legend()
        plt.show()
        #"""DENSITY"""
        plt.plot(mass_flow_range*1e3, density_g_equil, label = "density_g_equil")
        plt.plot(mass_flow_range*1e3, density_g_neq, label = "density_g_neq")
        plt.ylabel(r"$density_g$ [K]")
        plt.xlabel("mass flow [g/s]")
        plt.xlim(mass_flow_range[0]*1e3, mass_flow_range[-1]*1e3)
        plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
        plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
        # Make the minor ticks and gridlines show.
        plt.minorticks_on()
        plt.legend()
        plt.show()
        #"""VISCOSITY"""
        plt.plot(mass_flow_range*1e3, viscosity_g_equil, label = "viscosity_g_equil")
        plt.plot(mass_flow_range*1e3, viscosity_g_neq, label = "viscosity_g_neq")
        plt.ylabel(r"$viscosity_g$ [K]")
        plt.xlabel("mass flow [g/s]")
        plt.xlim(mass_flow_range[0]*1e3, mass_flow_range[-1]*1e3)
        plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
        plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
        # Make the minor ticks and gridlines show.
        plt.minorticks_on()
        plt.legend()
        plt.show()
        #"""MW"""
        plt.plot(mass_flow_range*1e3, MW_g_equil, label = "MW_g_equil")
        plt.plot(mass_flow_range*1e3, MW_g_neq, label = "MW_g_neq")
        plt.ylabel(r"$MW_g$ [K]")
        plt.xlabel("mass flow [g/s]")
        plt.xlim(mass_flow_range[0]*1e3, mass_flow_range[-1]*1e3)
        plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
        plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
        # Make the minor ticks and gridlines show.
        plt.minorticks_on()
        plt.legend()
        plt.show()


        def save_dq_study_to_excel(
            base_path,
            mass_flow_range,
            dH_range,
            T_g_equil,
            cp_g_equil,
            cv_g_equil,
            lambda_g_equil,
            density_g_equil,
            viscosity_g_equil,
            MW_g_equil,
            filename="dQ_study_thermal_properties_equil.xlsx"
        ):
            """
            Save dQ study equilibrium results into a single Excel sheet.
            """

            df = pd.DataFrame({
                "mass_flow_kg_s": mass_flow_range,
                "mass_flow_g_s": mass_flow_range * 1e3,
                "dH_removed_J_kg": dH_range,
                "T_g_equil_K": T_g_equil,
                "cp_g_equil_J_kgK": cp_g_equil,
                "cv_g_equil_J_kgK": cv_g_equil,
                "lambda_g_equil_W_mK": lambda_g_equil,
                "density_g_equil_kg_m3": density_g_equil,
                "viscosity_g_equil_Pa_s": viscosity_g_equil,
                "MW_g_equil_kg_kmol": MW_g_equil,
            })

            output_path = base_path / filename

            with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
                df.to_excel(writer, sheet_name="dQ_study_equil", index=False)

                ws = writer.sheets["dQ_study_equil"]
                ws.freeze_panes = "A2"

                for col in ws.columns:
                    max_length = 0
                    col_letter = col[0].column_letter
                    for cell in col:
                        try:
                            max_length = max(max_length, len(str(cell.value)))
                        except Exception:
                            pass
                    ws.column_dimensions[col_letter].width = max_length + 2

            print(f"Excel file saved to: {output_path}")    


        # def estimate_residence_time_mach(gas: ct.Solution,
        #                                 L_hx: float,
        #                                 mach: float = 0.3) -> tuple[float, float]:
        #     """
        #     Estimate gas residence time assuming a fixed Mach number.

        #     v_gas = M · a     where   a = sqrt(γ · r_specific · T)
        #     τ     = L_hx / v_gas

        #     Parameters
        #     ----------
        #     gas   : ct.Solution at HX inlet conditions (post-combustion state)
        #     L_hx  : heat exchanger length [m]
        #     mach  : flow Mach number (default 0.3)

        #     Returns
        #     -------
        #     tau   : residence time [s]
        #     v_gas : axial velocity [m/s]
        #     """
        #     gamma      = gas.cp / gas.cv                              # [-]
        #     r_specific = ct.gas_constant / gas.mean_molecular_weight  # J/kg/K  (MW in kg/kmol → /1000 not needed, Cantera uses SI)
        #     a          = np.sqrt(gamma * r_specific * gas.T)          # speed of sound [m/s]
        #     v_gas      = mach * a
        #     tau        = L_hx / v_gas
        #     return tau, v_gas

        # def finite_rate_hx_solve(gas_inlet: ct.Solution,
        #                         mdot: float,
        #                         removed_power: float,
        #                         residence_time: float,
        #                         n_steps: int = 500) -> ct.SolutionArray:
        #     """
        #     Finite-rate kinetics cooling of a combustion gas parcel in a heat exchanger.

        #     Plug-flow (Lagrangian parcel) model:
        #         IdealGasConstPressureReactor + Wall with prescribed heat flux.
        #     Heat is distributed uniformly along the HX length (constant q" assumption).

        #     Parameters
        #     ----------
        #     gas_inlet      : ct.Solution at post-combustion equilibrium state.
        #                     Pass a copy (gas.copy() or re-solve) — reactor modifies it.
        #     mdot           : hot gas mass flow [kg/s]
        #     removed_power  : total thermal power extracted over full HX [W]
        #     residence_time : gas parcel dwell time in HX [s]  →  τ = L_hx / v_gas
        #     n_steps        : integration steps (500–1000 for diesel surrogate)

        #     Returns
        #     -------
        #     states : ct.SolutionArray  (n_steps+1 rows)
        #         Access as states.T, states.cp, states.density, etc.
        #         Extra column: states.t  [s] (time = axial position proxy)
        #     """

        #     # ── Heat removal rate ────────────────────────────────────────────────────
        #     # Total specific enthalpy removal : dh = removed_power / mdot  [J/kg]
        #     # Distributed uniformly in time   : dh/dt = removed_power / (mdot·τ)  [W/kg]
        #     q_dot_per_kg = removed_power / (mdot * residence_time)   # W/kg

        #     # ── Reactor setup ────────────────────────────────────────────────────────
        #     r = ct.IdealGasConstPressureReactor(gas_inlet)   # constant-pressure PFR element

        #     # Cold sink — temperature is irrelevant because we prescribe flux, not U·ΔT
        #     sink_gas = ct.Solution(gas_inlet.source)
        #     sink_gas.TPX = 300.0, gas_inlet.P, gas_inlet.X
        #     sink = ct.Reservoir(sink_gas)

        #     # Wall: Q_dot [W] = heat_flux [W/m²] × area [m²]
        #     # Trick: set area = r.mass  →  heat_flux [W/m²] numerically equals q_dot_per_kg [W/kg]
        #     #        so  Q_dot_out = q_dot_per_kg × r.mass  ✓
        #     # r.mass is constant (no inlet/outlet in ConstPressureReactor).
        #     # Positive heat_flux ⟹ energy leaves reactor r (hot gas) ✓
        #     wall = ct.Wall(r, sink, A=r.mass)
        #     wall.heat_flux = q_dot_per_kg          # W/m²  (≡ W/kg with the area trick)

        #     net = ct.ReactorNet([r])
        #     net.rtol = 1e-9    # tight tolerances — diesel surrogate kinetics are stiff
        #     net.atol = 1e-15

        #     # ── Integration ──────────────────────────────────────────────────────────
        #     states = ct.SolutionArray(gas_inlet, extra=['t'])
        #     states.append(r.thermo.state, t=0.0)

        #     dt = residence_time / n_steps
        #     try:
        #         for i in range(1, n_steps + 1):
        #             net.advance(i * dt)
        #             states.append(r.thermo.state, t=i * dt)
        #     except ct.CanteraError as exc:
        #         print(f"[finite_rate_hx_solve] Integration failed at step {i}/{n_steps}: {exc}")
        #         print(f"  Partial result returned  (last T = {r.T:.1f} K)")

        #     return states


        # # ── Geometry inputs ──────────────────────────────────────────────────────────
        # L_hx  = 0.5   # m  — only free parameter now, A_cross no longer needed
        # MACH  = 0.3
        # # ─────────────────────────────────────────────────────────────────────────────
        # T_g_fr  = np.zeros(len(dH_range))   # finite rate results
        # cp_fr   = np.zeros(len(dH_range))
        # lambda_fr = np.zeros(len(dH_range))
        # density_fr = np.zeros(len(dH_range))
        # viscosity_fr = np.zeros(len(dH_range))
        # MW_fr   = np.zeros(len(dH_range))

        # for i, mdot in enumerate(mass_flow_range):

        #     # ── Solve combustion to get fresh inlet gas ───────────────────────────
        #     comb = combustion_gas_solve(
        #         fuel=fuel, oxidizer="O2",
        #         OF=O__F, p0=pcomb,
        #         T_g_init=T_ig, T_inj_LOX=T_inj,
        #         chem_mech_path=chem_mech_path,
        #         Hv_fuel=Hv_fuel, Y_fuel=Y_fuel)
        #     comb.solve()
        #     gas_inlet = comb.phase   # post-combustion equilibrium state

        #     # ── Estimate residence time from HX geometry ──────────────────────────
        #     tau, v_gas = estimate_residence_time_mach(gas_inlet, L_hx, mach=MACH)

        #     print(f"mdot={mdot*1e3:.1f} g/s | v_gas={v_gas:.1f} m/s | τ={tau*1e3:.2f} ms")

        #     # ── Finite-rate solve ─────────────────────────────────────────────────
        #     states = finite_rate_hx_solve(
        #         gas_inlet     = gas_inlet,
        #         mdot          = mdot,
        #         removed_power = removed_power,
        #         residence_time= tau,
        #         n_steps       = 500
        #     )

        #     # Final state (end of HX)
        #     T_g_fr[i]     = states.T[-1]
        #     cp_fr[i]      = states.cp[-1]
        #     lambda_fr[i]  = states.thermal_conductivity[-1]
        #     density_fr[i] = states.density[-1]
        #     viscosity_fr[i]= states.viscosity[-1]
        #     MW_fr[i]      = states.mean_molecular_weight[-1]

        # # ── Plot three-way comparison ─────────────────────────────────────────────────
        # fig, axes = plt.subplots(2, 3, figsize=(14, 8))
        # fig.suptitle("Equilibrium vs Frozen vs Finite-Rate — HX exit", fontsize=13)

        # props = [
        #     (T_g_equil,    T_g_neq,    T_g_fr,    r"$T_g$ [K]",         "Temperature"),
        #     (cp_g_equil,   cp_g_neq,   cp_fr,     r"$c_p$ [J/kg·K]",    "Heat capacity"),
        #     (lambda_g_equil,lambda_g_neq,lambda_fr,r"$\lambda$ [W/m·K]","Conductivity"),
        #     (density_g_equil,density_g_neq,density_fr,r"$\rho$ [kg/m³]","Density"),
        #     (viscosity_g_equil,viscosity_g_neq,viscosity_fr,r"$\mu$ [Pa·s]","Viscosity"),
        #     (MW_g_equil,   MW_g_neq,   MW_fr,     r"$MW$ [kg/kmol]",    "Mol. weight"),
        # ]

        # x = mass_flow_range * 1e3
        # for ax, (eq, frz, fr, ylabel, title) in zip(axes.flat, props):
        #     ax.plot(x, eq,  label="Equilibrium",  color="tab:blue")
        #     ax.plot(x, frz, label="Frozen",       color="tab:red",    linestyle="--")
        #     ax.plot(x, fr,  label="Finite rate",  color="tab:green",  linestyle="-.")
        #     ax.set_xlabel("ṁ [g/s]");  ax.set_ylabel(ylabel)
        #     ax.set_title(title);  ax.grid(True, alpha=0.4);  ax.legend(fontsize=8)

        # plt.tight_layout()
        # plt.show()




        save_dq_study_to_excel(
            base_path=base_path,
            mass_flow_range=mass_flow_range,
            dH_range=dH_range,
            T_g_equil=T_g_equil,
            cp_g_equil=cp_g_equil,
            cv_g_equil=cv_g_equil,
            lambda_g_equil=lambda_g_equil,
            density_g_equil=density_g_equil,
            viscosity_g_equil=viscosity_g_equil,
            MW_g_equil=MW_g_equil,
            filename="dQ_study_thermal_properties_equil.xlsx"
        )        #%%