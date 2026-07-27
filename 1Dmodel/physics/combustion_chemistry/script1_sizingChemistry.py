#%%
import cantera as ctr 
from pathlib import Path
import numpy as np 
import matplotlib.pyplot as plt
from pathlib import Path

ctr.suppress_thermo_warnings()

base_path = Path(__file__).parent


#%%
# Jet A (Kerosene) species "POSF10325"
fuel = "POSF10325"
JetA = "C:/Users/raubry/Desktop/Sizing codes/HX-combustor/A2highT.yaml"

chem_mech_path = base_path / "llnl_gasoline_323.yaml"
Y_fuel = {'C2H5OH': np.float64(0.04648289377205732), 'C6H12-1': np.float64(0.059442690615179886),
        'C5H10-1': np.float64(0.021229532362564234), 'IC4H8': np.float64(0.01132241726003426),
        'IC8H18': np.float64(0.5936889008300161), 'NC7H16': np.float64(0.2678335651601482)}
Hv_fuel = 313.664e3 # J/kg

# enthalpy of vaporization
Hv_JetA = 0.36e6 #J/kg  https://web.stanford.edu/group/haiwanglab/HyChem/approach/Report_Jet_Fuel_Thermochemical_Properties_v6.pdf
phase_chem_JetA = ctr.Solution(JetA) # instantiating the gas phase
phase_chem_JetA_2 = ctr.Solution(JetA) # instantiating the gas phase

# oxygen injection pressure and temperature
T_inj = 110 # K
p_inj = 1e5 # Pa, #! assume orifice dP

pcomb = 1e5

#%%

OF_list = np.linspace(0.01, 2.85, 100, endpoint=True)
T_g = np.zeros(len(OF_list))
cp_g = np.zeros(len(OF_list))
lambda_g = np.zeros(len(OF_list))
density_g = np.zeros(len(OF_list))

for i in range(len(OF_list)):

    O = OF_list[i]/(1 + OF_list[i])
    F = 1 - O

    T_mix = 1000 # 

    # evaporated fuel
    # Hv defined above


    """ENERGY LOSS LOX-->GOX"""

    # LOX TO GOX
    O2 = ctr.Oxygen() # import oxygen 
    O2.TP = T_inj, pcomb # set O2 to injection thermo
    h_ox_inj = O2.enthalpy_mass # fetch O2 onjection enthalpy
    O2.TP = T_mix, pcomb # set O2 to reaction state, assumed at T_pyro
    h_ox_pyro = O2.enthalpy_mass
    dH_ox = h_ox_pyro-h_ox_inj

    dH_tot = 0 #O*dH_ox + F*Hv_JetA


    """CHEMICAL EQUILIBRIUM CALCULATION"""


    phase_chem_JetA.TPY = T_mix, pcomb,  "O2" + ':' + str(OF_list[i]) + ',' +  fuel + ': 1' 

    phase_chem_JetA.equilibrate('HP') # finding chem equilibrium

    phase_chem_JetA.HP = phase_chem_JetA.enthalpy_mass-dH_tot, pcomb # remove energy loss
    phase_chem_JetA.equilibrate('HP') # re-equilibriate mixture

    T_g[i] = phase_chem_JetA.T
    cp_g[i] = phase_chem_JetA.cp
    lambda_g[i] = phase_chem_JetA.thermal_conductivity
    density_g[i] = phase_chem_JetA.density_mass

#%%

OF_max_Tg = OF_list[np.argmax(T_g)]
OF_st = 3.42

plt.plot(OF_list, T_g)
plt.ylabel(r"$T_g$ [K]")
plt.xlabel("O/F")
plt.xlim(OF_list[0], OF_list[-1])
plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
# Show the minor grid as well. Style it in very light gray as a thin,
# plt.axvline(x=OF_max_Tg, color='orange', linestyle='--', linewidth=1, label="OF@T_max")
# plt.axvline(x=OF_st, color='blue', linestyle='--', linewidth=1, label="OF@stoich")
# dotted line.
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
# Make the minor ticks and gridlines show.
plt.minorticks_on()
# plt.legend()
plt.show()

plt.plot(OF_list, cp_g/1e3)
plt.ylabel(r"$cp_g$ [kJ/kg-K]")
plt.xlabel("O/F")
plt.xlim(OF_list[0], OF_list[-1])
plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
# plt.axvline(x=OF_max_Tg, color='orange', linestyle='--', linewidth=1, label="OF@T_max")
# plt.axvline(x=OF_st, color='blue', linestyle='--', linewidth=1, label="OF@stoich")
# Show the minor grid as well. Style it in very light gray as a thin,
# dotted line.
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
# Make the minor ticks and gridlines show.
plt.minorticks_on()
# plt.legend()

plt.show()

plt.plot(OF_list, lambda_g)
plt.ylabel(r"$\lambda_g$ [W/m-K]")
plt.xlabel("O/F")
plt.xlim(OF_list[0], OF_list[-1])
plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
# plt.axvline(x=OF_max_Tg, color='orange', linestyle='--', linewidth=1, label="OF@T_max")
# plt.axvline(x=OF_st, color='blue', linestyle='--', linewidth=1, label="OF@stoich")

# Show the minor grid as well. Style it in very light gray as a thin,
# dotted line.
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
# Make the minor ticks and gridlines show.
plt.minorticks_on()
# plt.legend()

plt.show()

plt.plot(OF_list, density_g)
plt.ylabel(r"$\rho$ [kg/m3]")
plt.xlabel("O/F")
plt.xlim(OF_list[0], OF_list[-1])
plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
# plt.axvline(x=OF_max_Tg, color='orange', linestyle='--', linewidth=1, label="OF@T_max")
# plt.axvline(x=OF_st, color='blue', linestyle='--', linewidth=1, label="OF@stoich")

# Show the minor grid as well. Style it in very light gray as a thin,
# dotted line.
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
# Make the minor ticks and gridlines show.
plt.minorticks_on()
# plt.legend()

plt.show()

#%% mixing with Helium

perc_He_list = np.linspace(0, 0.9, 10)
OF_2_list = np.linspace(3, 40, 50)

density_2 = np.zeros((len(OF_2_list), len(perc_He_list)))
T_2 = np.zeros((len(OF_2_list), len(perc_He_list)))

pressurant  = "HE"
pcomb_2 = 70e5 

for i in range(len(perc_He_list)):
    for j in range(len(OF_2_list)):


        left_gas = 1-perc_He_list[i]

        O = OF_2_list[j]/(1 + OF_2_list[j])
        F = 1 - O

        perc_O = O * left_gas/(O + F)
        perc_F = F * left_gas/(O + F)


        T_mix = 800 # 

        # evaporated fuel
        # Hv defined above

        """ENERGY LOSS LOX-->GOX"""

        # LOX TO GOX
        O2 = ctr.Oxygen() # import oxygen 
        O2.TP = T_inj, pcomb_2 # set O2 to injection thermo
        h_ox_inj = O2.enthalpy_mass # fetch O2 onjection enthalpy
        O2.TP = T_mix, pcomb_2 # set O2 to reaction state, assumed at T_pyro
        h_ox_pyro = O2.enthalpy_mass
        dH_ox = h_ox_pyro-h_ox_inj

        dH_tot = perc_O*dH_ox + perc_F*Hv_JetA


        """CHEMICAL EQUILIBRIUM CALCULATION"""


        phase_chem_JetA_2.TPY = T_mix, pcomb_2,  "O2" + ':' + str(perc_O) + ',' +  fuel + ':' + str(perc_F) + ',' + pressurant + ':' + str(perc_He_list[i])
        phase_chem_JetA_2

        phase_chem_JetA_2.equilibrate('HP') # finding chem equilibrium

        phase_chem_JetA_2.HP = phase_chem_JetA_2.enthalpy_mass-dH_tot, pcomb_2 # remove energy loss
        phase_chem_JetA_2.equilibrate('HP') # re-equilibriate mixture

        # """ Remove Helium heating energy"""
        # dH_he = perc_He_list[i] * 5200*(phase_chem_JetA_2.T - 100)
        # phase_chem_JetA_2.HP = phase_chem_JetA_2.enthalpy_mass-dH_he, pcomb_2 # remove energy loss
        # phase_chem_JetA_2.equilibrate('HP') # re-equilibriate mixture
        

        density_2[j, i] = phase_chem_JetA_2.density
        T_2[j, i] = phase_chem_JetA_2.T

#%%

for k in range(len(perc_He_list)):

    plt.plot(OF_2_list, density_2[:,k], label = f"%He={np.round(perc_He_list[k],2)}")

plt.ylabel(r"$\rho$ [kg/m3]")
plt.xlabel("O/F")
plt.xlim(OF_2_list[0], OF_2_list[-1])
plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
# Show the minor grid as well. Style it in very light gray as a thin,
# dotted line.
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
# Make the minor ticks and gridlines show.
plt.minorticks_on()
# plt.legend()
plt.show()


for k in range(len(perc_He_list)):

    plt.plot(OF_2_list, T_2[:,k], label = f"%He={np.round(perc_He_list[k],2)}")

plt.ylabel(r"$T_g$ [K]")
plt.xlabel("O/F")
plt.xlim(OF_2_list[0], OF_2_list[-1])
plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
# Show the minor grid as well. Style it in very light gray as a thin,
# dotted line.
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
# Make the minor ticks and gridlines show.
plt.minorticks_on()
# plt.legend()
plt.show()


#%%