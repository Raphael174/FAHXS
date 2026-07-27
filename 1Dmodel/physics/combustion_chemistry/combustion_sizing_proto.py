""" 
@ author : Raphaël Aubry

Script for the pre-sizing of the combustor

-compute characteristic area of torch ignitor to reach or avoid sonic throat
"""
#%%
import numpy as np 
from combustion_gas import choose_fuel, combustion_gas_solve
from scipy.constants import gas_constant, g
from CoolProp.CoolProp import PropsSI
#%% input data

mass_flow_tot = 3e-3 # kg/s, approximate total, with +/- 20% margin
OF_tot = 1 #3.682

fuel = "gasoline-E10"
density_fuel = 750 #kg/m3, avg at 20°C, for SP95
oxidizer = "O2"
massTorch__masstot = 0.0
mass_flow_torch = massTorch__masstot*mass_flow_tot
OF_torch = 0.3
Dh_torch = 50e-3
Dt_torch = 10e-3
At_torch = np.pi*Dt_torch**2/4

spark_temp = 500 + 273 # K 
T_inj_O2 = -20+273 # K


mass_flow_main = mass_flow_tot - mass_flow_torch
Dt_cc = 20e-3
At_cc = np.pi*Dt_cc**2/4

mass_flow_ox_main = mass_flow_main*OF_tot/(1 + OF_tot)*1.1
mass_flow_fu_main = mass_flow_main/(1 + OF_tot)

p0_burner = 2 # bar, total pressure of burner
p0_torch = 2

# pipe limit velocities
Vmax_main_pipe_O2 = 10 # m/s, at 40 bars 
Vmax_secondary_O2 = 20 # m/s, at 30 bar

# injector sizing
p_line =6
N_main = 3 # number of injector ports in main chamber, both fuel and oxidizer
mass_flow_ox_main_unit = mass_flow_ox_main/N_main
mass_flow_fu_main_unit = mass_flow_fu_main/N_main

#gasoline injector cc/min rating
Q_fu_ccMin = mass_flow_fu_main_unit/density_fuel * 100**3 * 60

# assumed O2 orifice discharge coefficient
Cd = 0.6

cp_O2 = PropsSI('CPMASS','T', T_inj_O2,'P',p_line*1e5 ,oxidizer)
cv_O2 = PropsSI('CVMASS','T', T_inj_O2,'P',p_line*1e5 ,oxidizer)
gamma_O2 = cp_O2/cv_O2
density_O2 = PropsSI('D','T', T_inj_O2,'P',p_line*1e5 ,oxidizer)

# estimating O2 injection velocity at burner pressure
Dh_O2_exit = 12e-3 # sortie en G1/4"
A_O2_exit = np.pi*Dh_O2_exit**2/4
density_O2_inj = PropsSI('D','T', T_inj_O2,'P',p0_burner*1e5 ,oxidizer)
V_O2_inj = mass_flow_ox_main_unit/(density_O2_inj*A_O2_exit)

# O2 min pipe diameters according to max velocities required 
density_O2_main_maxP = PropsSI('D','T', T_inj_O2,'P',40e5 ,oxidizer)
density_O2_sec_maxP = PropsSI('D','T', T_inj_O2,'P',20e5 ,oxidizer)

Dh_pipe_O2_main = np.sqrt(4/np.pi * mass_flow_ox_main/(density_O2_main_maxP*Vmax_main_pipe_O2))
Dh_pipe_O2_sec = np.sqrt(4/np.pi * mass_flow_ox_main/N_main/(density_O2_sec_maxP*Vmax_secondary_O2))
#%%

# combustion chamber
def Area_sonic ( m_dot, gamma, pc, R, Tc):
    Gam = np.sqrt(gamma)*(2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
    return m_dot*np.sqrt(R*Tc)/(Gam*pc)

def p0_chamber(m_dot, gamma, At, R, Tc):
    Gam = np.sqrt(gamma)*(2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
    return m_dot*np.sqrt(R*Tc)/(Gam*At)

def exit_velocity (pe, p0, gamma, R, Tc):
    return np.sqrt(2*gamma/(gamma-1)*R*Tc*(1 - (pe/p0)**((gamma-1)/gamma)))

def pressure_at_Mach (p0, gamma, Ma):
    return p0/(1 + (gamma-1)/2*Ma**2)**(gamma/(gamma-1))

# injectors

# choked orifice
def diameter_orifice_choked (Cd, m_dot, gamma, p_upstream, rho_upstream):
    """ 
    Alessandro de Iaco Veris - Fundamental Concepts of Liquid-Propellant Rocket Engines (2021, Springer)
    p516
    """
    term1 = np.sqrt(gamma*p_upstream*rho_upstream * (2/(gamma+1))**((gamma+1)/(gamma-1)))
    return np.sqrt(4*m_dot/(Cd*np.pi*term1))

def critical_pressure_ratio (gamma):
    """
    critical downstream pressure/upstream total pressure
    Ideal gas assumption
    """
    return (2/(gamma+1))**(gamma/(gamma-1))

#%% injector 
print(" ")
print("INJECTION")
print(f"mass_flow_tot ={mass_flow_tot*1e3}g/s | %torch_mass_flow ={massTorch__masstot}%")
print(f"mass_flow_fu_main={mass_flow_fu_main*1e3} ")
print("")
print(f"Diameters min size, main = {Dh_pipe_O2_main*1e3} | sec = {Dh_pipe_O2_sec*1e3} mm")
print(" ")
print(f"Gasoline injector rating = {Q_fu_ccMin} cc/min")
print(" ")
print(f"Injector elements = {N_main}")
print(f"mass_flow_ox_unit={mass_flow_ox_main_unit*1e3}g/s | mass_flow_fu_unit={mass_flow_fu_main_unit*1e3}g/s")
critical_pressure_ratio_O2 = critical_pressure_ratio(gamma=gamma_O2)
p_upstream_min_choke = p0_burner/critical_pressure_ratio_O2
d_orifice = diameter_orifice_choked(Cd=Cd, m_dot=mass_flow_ox_main_unit, gamma=gamma_O2, p_upstream=p_line*1e5, rho_upstream=density_O2)
p_down = critical_pressure_ratio(gamma_O2)*p_line
print(" ")
print("Assuming orifice downstream pressure = combustor target pressure")
print(f"Minimum upstream pressure for O2 choked orifice = {p_upstream_min_choke} bars")

print(f"Orifice diameter based on unitary O2 flow m_ox/{N_main}-m_ox_ig and upsream line pressure of {p_line} bar = {d_orifice*1e3}mm, p_downstream={p_down} bar")

print(f"velocity O2 into chamber = {V_O2_inj}m/s")

#%% Combustion chamber
"""
TORCH
"""
chem_mech_path, Y_fuel, Hv_fuel = choose_fuel(fuel)

combustion_object_torch = combustion_gas_solve(
                fuel=fuel, oxidizer=oxidizer,
                OF=OF_torch,
                p0=p0_burner*1e5,
                T_g_init=spark_temp,
                T_inj_LOX=T_inj_O2,
                chem_mech_path=chem_mech_path, Hv_fuel=Hv_fuel, Y_fuel=Y_fuel)
combustion_object_torch.solve()
combustion_object_torch.phase.T

T_torch = combustion_object_torch.phase.T
R_torch = gas_constant*1e3/combustion_object_torch.phase.mean_molecular_weight
gamma_torch = combustion_object_torch.phase.cp/combustion_object_torch.phase.cv 

A_sonic_torch = Area_sonic(m_dot=massTorch__masstot*mass_flow_tot, gamma=gamma_torch,
                           pc=p0_torch*1e5, R=R_torch, Tc=T_torch)
D_sonic_torch = np.sqrt(4*A_sonic_torch/np.pi)

p0_at_At = p0_chamber(m_dot=massTorch__masstot*mass_flow_tot, gamma=gamma_torch,
                           At=At_torch, R=R_torch, Tc=T_torch)

pe_torch = pressure_at_Mach(p0=p0_torch, gamma=gamma_torch, Ma=1)
Ve_torch = exit_velocity(pe=pe_torch, p0=p0_torch, gamma=gamma_torch, R=R_torch, Tc=T_torch)

print(" ")
print("IGNITION")
print(f"D_sonic_torch = {D_sonic_torch*1e3} mm @ p0_burner={p0_burner}bar")
print(f"T_torch = {T_torch-273.15} C°")
print(f"gamma_torch = {gamma_torch}")
print(f"m_dot_torch = {massTorch__masstot*mass_flow_tot*1e3} g/s")
print(f"R_torch = {R_torch} ")
print(f"O/F_torch = {OF_torch}")

print(f"pe={pe_torch/1e5} bar, @ p0={p0_torch} bar")
print(f"Ve={Ve_torch} m/s, @ p0={p0_torch} bar")

""" 
MAIN CC
"""

combustion_object_cc = combustion_gas_solve(
                fuel=fuel, oxidizer=oxidizer,
                OF=OF_tot,
                p0=p0_burner*1e5,
                T_g_init=spark_temp,
                T_inj_LOX=T_inj_O2,
                chem_mech_path=chem_mech_path, Hv_fuel=Hv_fuel, Y_fuel=Y_fuel)
combustion_object_cc.solve()
combustion_object_cc.phase.T

T_cc = combustion_object_cc.phase.T
R_cc = gas_constant*1e3/combustion_object_cc.phase.mean_molecular_weight
gamma_cc = combustion_object_cc.phase.cp/combustion_object_cc.phase.cv 

A_sonic_cc = Area_sonic(m_dot=mass_flow_main, gamma=gamma_cc,
                           pc=p0_burner*1e5, R=R_cc, Tc=T_cc)
D_sonic_cc = np.sqrt(4*A_sonic_cc/np.pi)

p0_at_At = p0_chamber(m_dot=mass_flow_main, gamma=gamma_cc,
                           At=At_cc, R=R_cc, Tc=T_cc)

pe_cc = pressure_at_Mach(p0=p0_at_At, gamma=gamma_cc, Ma=1)
Ve_cc = exit_velocity(pe=pe_cc, p0=p0_at_At, gamma=gamma_cc, R=R_cc, Tc=T_cc)

print(" ")
print("CHAMBER")
print(f"D_sonic_cc = {D_sonic_cc*1e3} mm")
print(f"T_cc = {T_cc-273.15} C°")
print(f"gamma_cc = {gamma_cc}")
print(f"m_dot_main = {mass_flow_main*1e3} g/s")
print(f"R_cc = {R_cc} ")
print(f"O/F_torch = {OF_tot}")

print(f"@ Dt = {Dt_cc*1e3} mm")
print(f"p0={p0_at_At/1e5} | pe={pe_cc/1e5} bar")
print(f"Ve={Ve_cc} m/s, @ p0={p0_burner} bar")


Thrust = mass_flow_tot*Ve_cc + (pe_cc-1e5)*At_cc
thrust__g = Thrust/g
print(f"Thrust = {Thrust} N")
print(f"Thrust mass = {thrust__g} kg")

#%%