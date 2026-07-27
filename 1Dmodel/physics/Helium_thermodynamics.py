"""
@ author : Raphaël Aubry

Standalone exploratory script (Jupyter-style #%% cells) plotting He
compressibility factor Z vs T,p. Not imported by any solver or test — audited
2026-07-13 during the liquid-coolant physics restructure and confirmed to have
no importers; left in place rather than moved into physics/gas_flow/, since it
is offline analysis, not a reusable module.
"""
#%%
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_inline.backend_inline import set_matplotlib_formats
set_matplotlib_formats('svg')


#%%
fluid = "HELIUM"

#%%
range_T = np.array([50, 200, 400, 600, 800])
range_p = np.linspace(70, 90, 50, endpoint=True)*1e5

cp = np.zeros((len(range_p),len(range_T)))
k = np.zeros((len(range_p),len(range_T)))
mu = np.zeros((len(range_p),len(range_T)))
rho = np.zeros((len(range_p),len(range_T)))
cv = np.zeros((len(range_p),len(range_T)))
Z = np.zeros((len(range_p),len(range_T)))

for i in range(len(range_T)):
    for j in range(len(range_p)): 
        cp[j,i] = PropsSI('CPMASS','T', range_T[i],'P',range_p[j],fluid)
        cv[j,i] = PropsSI('CVMASS','T', range_T[i],'P',range_p[j],fluid)
        k[j,i]  = PropsSI('L','T',range_T[i],'P',range_p[j],fluid)
        mu[j,i] = PropsSI('V','T',range_T[i],'P',range_p[j],fluid)
        Z[j,i]  = PropsSI('Z','T',range_T[i],'P',range_p[j],fluid)

#%%

for i in range(len(range_T)):
    plt.plot(range_p/1e5, Z[:,i], label=f"{range_T[i]}K")

plt.xlabel(r"$p$ [bar]")
plt.ylabel(r"$p/\rho R T$")
plt.xlim(range_p[0]/1e5, range_p[-1]/1e5)
plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
plt.grid(which='minor', color='#EEEEEE', linestyle='--', linewidth=0.8)
plt.minorticks_on()
plt.legend()
plt.plot()

#%%

range_T_ = np.linspace(50, 800, 50, endpoint=True)
range_p_ = np.array([70, 75, 80, 85, 90])*1e5

cp = np.zeros((len(range_T_),len(range_p_)))
k = np.zeros((len(range_T_),len(range_p_)))
mu = np.zeros((len(range_T_),len(range_p_)))
rho = np.zeros((len(range_T_),len(range_p_)))
cv = np.zeros((len(range_T_),len(range_p_)))
Z = np.zeros((len(range_T_),len(range_p_)))

for i in range(len(range_p_)):
    for k in range(len(range_T_)): 
        cp[k,i] = PropsSI('CPMASS','T', range_T_[k],'P',range_p_[i],fluid)
        cv[k,i] = PropsSI('CVMASS','T', range_T_[k],'P',range_p_[i],fluid)
        #k[k,i]  = PropsSI('L','T',range_T_[k],'P',range_p_[i],fluid)
        mu[k,i] = PropsSI('V','T',range_T_[k],'P',range_p_[i],fluid)
        Z[k,i]  = PropsSI('Z','T',range_T_[k],'P',range_p_[i],fluid)

#%%

for i in range(len(range_p_)):
    plt.plot(range_T_, Z[:,i], label=f"{np.round(range_p_[i]/1e5, 0)} bar")

plt.xlabel(r"$T$ [K]")
plt.ylabel(r"$p/\rho R T$")
plt.xlim(range_T_[0], range_T_[-1])
plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
plt.grid(which='minor', color='#EEEEEE', linestyle='--', linewidth=0.8)
plt.minorticks_on()
plt.legend()
plt.plot()

#%%