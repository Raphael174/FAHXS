"""
@author : Raphael Aubry

Here we take the scripts for liquid and gas orifice sizing and combien them to 
reach an approriate momentum ratio 

Also we compute the estimated droplet size of doublets w/o gas flow, which should
underestimate spray droplet size 
"""

import matplotlib.pyplot as plt
import numpy as np 
from CoolProp.CoolProp import PropsSI
from scipy.constants import gas_constant
from liquid_impinging_orifice_main import straight_orifice_flow_LIQUID
from straight_gas_orifice_main import straight_orifice_flow_GAS
from impinging_injection_correlations import Webber_jet, jet_breakup_length_single_jet_GrantMiddleman, d10_like_doublet_Sweenie, Webber_transition_breakup_mode_like_doublet_Sweenie, sheet_breakup_length_like_doublet_Sweenie, length_impingement


#* Fuel data 
# https://iopscience.iop.org/article/10.1088/1742-6596/1473/1/012042/pdf

#! Chamber operating conditions
p_chamber = 1e5 # max 3 bar
mass_flow_tot = 50e-3 # kg/s, approximate total, with +/- 20% margin
OF_tot = 3.682
oxidizer = "O2"
fuel = "gasoline-E10"

#! Ignitor operating conditions
massTorch__masstot = 0.0 # % total propellant through torch
mass_flow_torch = massTorch__masstot*mass_flow_tot
OF_torch = 0.3

#! Mass share chamber and ignitor 
mass_flow_main = mass_flow_tot - mass_flow_torch

#!-------------------------
#* MOMENTUM RATIO GAS/LIQUID
rhogUg__rhoLUL = 1.5
#*
#!-------------------------

#! Injector (LIQUID) operating conditions
dP_inj_liquid= 1.5e5 #p_chamber*dP__p0
T_fuel = 18 + 273
density_fuel = 730 #kg/m3, avg at 20°C, for SP95
kinematic_viscosity_fuel = 0.5383*1e-6
mass_flow_fu_main = mass_flow_main/(1 + OF_tot)
surface_tension_liquid = 0.018 # N-m

jets__element=2
N_elements = 3
l_o_L = 3.45e-3
entrance_type="tilted" # "tilted", "chamfer", "round"
alpha_entrance = 60 # jet angle
tolerance=1e-6

#* cavitation
MoS_cav = 1.
#y ≈ a * exp(b * (x - x0)) 
#Fit params: {'a': 80.52932489462026, 'b': 0.02933396938046209, 'x0': 337.63151999999997, 'r2': 0.9501150084945873}

def vapour_pressure_gasoline(T):
    """ 
    Approximate curve for gasoline vapour pressure
    https://sci-hub.se/10.1016/j.enconman.2015.10.081

    Returns pressure in Pa with input Kelvins
    """
    params = {'a': 80.52932489462026, 'b': 0.02933396938046209, 'x0': 337.63151999999997, 'r2': 0.9501150084945873}
    return 1000*params['a'] * np.exp(params['b'] * (T - params['x0'])) 

#! Injector (GAS) operating conditions
mass_flow_ox_main = mass_flow_main*OF_tot/(1 + OF_tot)
dP_inj_gas= 5e5 #2e5 #p_chamber*dP__p0 #* keep low for Ujet<400m/s
p_plenum=p_chamber+dP_inj_gas
T_plenum_O2 = 15 + 273.15 #* quite an assumption, but low impact on results
cp_O2 = PropsSI('CPMASS','T', T_plenum_O2,'P',p_plenum ,oxidizer)
cv_O2 = PropsSI('CVMASS','T', T_plenum_O2,'P',p_plenum ,oxidizer)
viscosity_O2 = PropsSI('V','T', T_plenum_O2,'P',p_plenum ,oxidizer)
density_O2 = PropsSI('D','T', T_plenum_O2,'P',p_plenum ,oxidizer)
R_O2 = 1e3*gas_constant/31.999

radial_distance_to_center = 50e-3
entrance_type = "chamfer"
beta_entrance_gas = 40
l_in_gas = 5e-3
l_o_G = 10e-3


#* Computations for liquid jets
liquid_orifice = straight_orifice_flow_LIQUID(  dP_inj=dP_inj_liquid, p_chamber=p_chamber, mass_flow_total=mass_flow_fu_main, 
                                            N_elements=N_elements, l_o=l_o_L, jets__element=jets__element, alpha_entrance=alpha_entrance, entrance_type=entrance_type)
liquid_orifice.solve()

print("LIQUID JETS SIZING : ")
print(f"iterations={liquid_orifice.i}")
print(f"total main fuel flow ={mass_flow_fu_main}")
print(f"C_d={liquid_orifice.C_d} | xi={liquid_orifice.xi}, xi_vortex={liquid_orifice.xi_1_c}, xi_int={liquid_orifice.x_in}, xi_fric={liquid_orifice.xi_fric}")
print(f"alpha={alpha_entrance}")
print(f"U_1={liquid_orifice.U_1} | U_2={liquid_orifice.U_2} m/s")
print(f"dP={dP_inj_liquid/1e5} bar")
print(f"do={liquid_orifice.d_o*1e3} | lo={liquid_orifice.l_o*1e3} mm")
print(f"Reynolds entrance={liquid_orifice.Re}")
print(f"eta_cav={liquid_orifice.eta}")
print(f"p_n = p_s - dP/eta = {liquid_orifice.p_n} >= {vapour_pressure_gasoline(T=T_fuel)}")

#cavitation computation
p0_needed = (vapour_pressure_gasoline(T=T_fuel) + dP_inj_liquid/liquid_orifice.eta)*MoS_cav/1e5
print(f"po needed to avoid cavitation={p0_needed} bar")

print(" ")
print("gas-liquid jets")
print(f"rho_L*U_2^2={liquid_orifice.momentum} kg/m-s²")
print(" ")
error_liquid = 0.015
liquid_orifice.solve_mass_flow(do=liquid_orifice.d_o+error_liquid*1e-3)
print(f"+{error_liquid}mm error on do--> mass={liquid_orifice.mass_flow_inj_*1e3} g/s v.s. mass={liquid_orifice.mass_flow_inj*1e3}")


# jet momentum 


#* Computations for gas orifice
gas_orifice = straight_orifice_flow_GAS(  dP_inj=dP_inj_gas, p_chamber=p_chamber, mass_flow_total=mass_flow_ox_main, 
                                            N_elements=N_elements, entrance_type=entrance_type, beta_entrance=beta_entrance_gas,
                                            gas_density_plenum=density_O2, gas_R=R_O2, gas_gamma=cp_O2/cv_O2, gas_T_plenum=T_plenum_O2, gas_dynamic_viscosity=viscosity_O2,
                                            l_in=l_in_gas, l_o=l_o_G, distance_to_center=radial_distance_to_center)
gas_orifice.solve()

print("---------------------")
print("GAS JET SIZING")
print(f"iterations={gas_orifice.i}")
print(f"total main oxygen flow ={mass_flow_ox_main}")
print(f"C_d={gas_orifice.C_d} | xi={gas_orifice.xi}, xi_vortex=assumed {0}, xi_int={gas_orifice.xi_in}, xi_fric={gas_orifice.xi_fric}")
print(f"U_1={gas_orifice.U_1} | U_2={gas_orifice.U_2} m/s")
print(f"dP={dP_inj_gas/1e5} bar")
print(f"do={gas_orifice.d_o*1e3} | lo={gas_orifice.l_o*1e3} mm")
print(f"Reynolds entrance={gas_orifice.Re}")
print(" ")
print("gas-liquid jets")
print(f"rho_L*U_2^2={gas_orifice.momentum} kg/m-s²")
# jet momentum 

#* Combined jets

jet_webber = Webber_jet(density_L=density_fuel, velocity_jet=liquid_orifice.U_2, l=liquid_orifice.d_o, surface_tension_L=surface_tension_liquid)
jet_breakup_length = jet_breakup_length_single_jet_GrantMiddleman(Weber_jet=jet_webber, do=liquid_orifice.d_o)


length_impingement = length_impingement(alpha=alpha_entrance, distance_jet_orifice=gas_orifice.d_o)
Webber_trans = Webber_transition_breakup_mode_like_doublet_Sweenie(theta=alpha_entrance, l_b=jet_breakup_length, l_i=length_impingement)
sheet_breakup_length = sheet_breakup_length_like_doublet_Sweenie(Weber_jet=jet_webber, Weber_transition=Webber_trans, do=liquid_orifice.d_o, theta=alpha_entrance, l_b=jet_breakup_length, l_i=length_impingement)
mean_spray_droplet_size = d10_like_doublet_Sweenie(Weber_jet=jet_webber, theta=alpha_entrance, l_b=jet_breakup_length, l_i=length_impingement)

print("---------------------")
print("COMBINED JETS")
print(f"GAS/LIQUID momentum={gas_orifice.momentum/liquid_orifice.momentum}")
print(f"length_impingement={length_impingement*1e3} mm")
print(f"break up length={jet_breakup_length*1e3} mm")
print(f"sheet_breakup_length={sheet_breakup_length*1e3} mm")
print(f"mean_spray_droplet_size={mean_spray_droplet_size*1e6} µm")



#%% study on the orifice diameters

study_1=True
if study_1==True:
    """ 
    Some things have been fixed and understood :
        -Chamber pressure @ 1atm as much as possible (low velocity exhaust better)
        -dP liquid maintained <= 1bar to avoid cavitation
        -dP gas ~ 2*dP_liquid for good momentum ratio gas/liquid
        -orifice length gas and liquid fixed throughout campaign, lo_L=3.45mm, lo_G=10mm (low impact on dP variation)
        -N_elements=3 good size for flow/combustion symmetry and for orifice size range
        -fluid temperatures assumed between 15-20C, sat pressure sized at 20C for fuel (oversize)
        -impinging jet droplets estimated all cases < 200 µm, so with O2 flow improved atomization expected
        -

        Now the goal is to find the orifice diameter ranges for gas and li
    """
    mass_flow_total = np.linspace(25, 150, 6, endpoint=True)*1e-3
    mass_fu = np.zeros(len(mass_flow_total))
    mass_ox = np.zeros(len(mass_flow_total))
    check_do_gas = np.zeros(len(mass_flow_total))
    check_do_liq = np.zeros(len(mass_flow_total))
    momentum_ratio = np.zeros(len(mass_flow_total))

    Cd_O2 = np.zeros(len(mass_flow_total))
    Cd_fu = np.zeros(len(mass_flow_total))

    OF_tot = 3.682
    massTorch__masstot = 0.0 # % total propellant through torch
    OF_torch = 0.3



    for i in range(len(mass_flow_total)):

        m_tot = mass_flow_total[i]
        mass_flow_torch = massTorch__masstot*m_tot
        mass_flow_main = m_tot - mass_flow_torch
        mass_flow_fu_main = mass_flow_main/(1 + OF_tot)
        mass_flow_ox_main = mass_flow_main*OF_tot/(1 + OF_tot)

        liquid_orifice_ = straight_orifice_flow_LIQUID(  dP_inj=dP_inj_liquid, p_chamber=p_chamber, mass_flow_total=mass_flow_fu_main, 
                                                    N_elements=N_elements, l_o=l_o_L, jets__element=jets__element, alpha_entrance=alpha_entrance, entrance_type=entrance_type)
        liquid_orifice_.solve()
        gas_orifice_ = straight_orifice_flow_GAS(  dP_inj=dP_inj_gas, p_chamber=p_chamber, mass_flow_total=mass_flow_ox_main, 
                                                    N_elements=N_elements, entrance_type=entrance_type, beta_entrance=beta_entrance_gas,
                                                    gas_density_plenum=density_O2, gas_R=R_O2, gas_gamma=cp_O2/cv_O2, gas_T_plenum=T_plenum_O2, gas_dynamic_viscosity=viscosity_O2,
                                                    l_in=l_in_gas, l_o=l_o_G, distance_to_center=radial_distance_to_center)
        gas_orifice_.solve()

        check_do_gas[i] = gas_orifice_.d_o*1e3
        check_do_liq[i] = liquid_orifice_.d_o*1e3
        mass_fu[i] = mass_flow_fu_main
        mass_ox[i] = mass_flow_ox_main
        momentum_ratio[i] = gas_orifice_.momentum/liquid_orifice_.momentum
        
        Cd_O2[i] = gas_orifice.C_d
        Cd_fu[i] = liquid_orifice.C_d_

    print("mass flows =", mass_flow_total)
    print("mass ox = ", mass_ox)
    print("mass fuel = ", mass_fu)
    print(("------------------------------------"))
    print("do_liquid =", check_do_liq)
    print("do_gas =", check_do_gas)
    print("------------------------------------")
    print("Cd_O2 = ", Cd_O2)
    print("Cd_fu = ", Cd_fu)
    print("------------------------------------")
    print("density O2 = ", density_O2, " density fuel = ", density_fuel, " | kg/m3")




    fig, ax1 = plt.subplots()

    # First axis: Res_g, Res_c, Res_w
    ax1.plot(mass_flow_total*1e3, check_do_gas, linestyle='-',label="do_ox", color="blue")
    ax1.plot(mass_flow_total*1e3, check_do_liq, linestyle='-', label="do_fu", color="orange")
    ax1.set_xlabel(r"$m_{tot}$ [g/s]")
    ax1.set_ylabel(r"$d_o$ [mm]")
    ax1.set_xlim(mass_flow_total[0]*1e3, mass_flow_total[-1]*1e3)
    ax1.ticklabel_format(useOffset=False) 
    ax1.grid(which='major', color='#DDDDDD', linewidth=0.8)
    ax1.grid(which='minor', color='#EEEEEE', linestyle='--', linewidth=0.8)
    ax1.minorticks_on()

    # Second axis: UA
    ax2 = ax1.twinx()
    ax2.plot(mass_flow_total*1e3, mass_fu*1e3, linestyle='--', label="m_fu", color="orange")
    ax2.plot(mass_flow_total*1e3, mass_ox*1e3, linestyle='--', label="m_ox", color="blue")
    ax2.set_ylabel(r"$m$ [g/s]")

    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')

    fig.show()
    plt.show()

    plt.plot(mass_flow_total*1e3, momentum_ratio)
    plt.xlim(mass_flow_total[0]*1e3, mass_flow_total[-1]*1e3)
    plt.ylabel(r"$\rho_g U_g^2 / \rho_l U_l^2$")
    plt.show()