"""
Liquid monopropellant jet-injector discharge coefficient per Yang (2005), Ch. 2
(II. Theory and Design of Liquid Monopropellant Jet Injectors).

Simple version: no 'typing' or '__future__' imports; clear variables and math only.
Works with the companion 'yang2005_digitized_curves.py' for ξ-inlet and ξ_1→c.

Definition:
    m_dot = C_d * A_n * sqrt(2 * rho * Δp) , with  C_d = ε / sqrt(1 + ξ_total)
and
    ξ_total = ξ_(1→c) + ξ_(c→2) + ξ_fr + ξ_m ,  ξ_fr = λ * (L/D) , λ ≈ 0.3164 Re^{-0.25}

Two regimes for ε:
    L/D ≥ 1.5  → ε ≈ 1.0 (reattached jet)
    L/D < 1.5  → user must provide ε (or allow a default), since free-jet forms
"""


#%%

from CoolProp.CoolProp import PropsSI
from straight_orifice_loss_factors import xi1c_from_fig22, xiin_from_fig23, xiin_from_fig24, xiin_from_fig25
import numpy as np 

#%%

# -----------------------------
# Core correlations
# -----------------------------

#! KEEP
def blasius_lambda(Re):
    """Darcy friction factor for smooth turbulent pipe; laminar fallback."""
    if Re <= 0:
        raise ValueError("Re must be positive.")
    elif Re < 2300.0:
        return 64.0 / Re
    else:
        return 0.3164 * Re ** (-0.25)


def choose_entrance_loss(entrance_type, beta, tau__d_o, alpha, l_in__d_o):
    if entrance_type == "tilted" :
        return xiin_from_fig25(alpha)
    elif entrance_type == "chamfer" :
        return xiin_from_fig23(beta_deg=beta, lin_over_di=l_in__d_o)
    elif entrance_type == "round" :
        return xiin_from_fig24(tau_over_di=tau__d_o)
    else:
        raise ValueError("unknown entrance type")


def lambda_2(gamma, p_chamber, p_plenum):
    """ 
        velocity coefficient lambda_2 for gas flow orifice according to p52 of 
        Method proposed by V. BAZAROV, V. YANG, P. PURI, in
        Liquid Rocket Thrust Chambers - V. Yang, M. Habiballah, M. Popp, J. Hulka - 2005
    """
    return np.sqrt((gamma+1)/(gamma-1) * ( 1 - (p_chamber/p_plenum)**((gamma-1)/gamma) ))

def Cstar(gamma, R, T_plenum):
    return np.sqrt(gamma*R*T_plenum)/ (gamma*np.sqrt( (2/(gamma+1))**((gamma+1)/(gamma-1) ) ) )

def q_lambda2(gamma, lambda_2):
    """ 
    parameter on page 52 of 
        Method proposed by V. BAZAROV, V. YANG, P. PURI, in
        Liquid Rocket Thrust Chambers - V. Yang, M. Habiballah, M. Popp, J. Hulka - 2005
    """
    return ((gamma+1)/2)**(1/(gamma-1)) * lambda_2*(1 + (gamma-1)/(gamma+1)**lambda_2**2)**(1/(gamma-1))

def gas_exit_velocity(gamma, R, T_plenum, p_chamber, p_plenum):
    """ 
    Isentropic gas exit velocity 
    """
    return np.sqrt(2*gamma/(gamma-1) *R*T_plenum*(1 - (p_chamber/p_plenum)**((gamma-1)/gamma)))

class straight_orifice_flow_GAS ():
    def __init__(self, 
                 dP_inj, mass_flow_total, N_elements, p_chamber,
                 #SP95
                 gas_density_plenum, gas_R, gas_gamma, gas_T_plenum, gas_dynamic_viscosity,
                 # 
                 d_o=7e-3, l_o=20e-3, l_in=5e-3, distance_to_center=50e-3,
                 beta_entrance =40, tau_entrance=1e-3,
                 entrance_type="tilted", # "tilted", "chamfer", "round"
                 tolerance=1e-6):
        """
        p51-53
        Method proposed by V. BAZAROV, V. YANG, P. PURI, in
        Liquid Rocket Thrust Chambers - V. Yang, M. Habiballah, M. Popp, J. Hulka - 2005

        Based on injector pressure drop, total mass flow required, 
        number of injector elements, fluid properties, compute : 
            orifice diameter
            exit jet velocity
            coefficient of discharge
        """
        #* expected working conditions
        self.dP_inj = dP_inj
        self.p_chamber = p_chamber
        self.mass_flow_total = mass_flow_total
        self.N_elements = N_elements #! number of impinging elements
        self.distance_to_center = distance_to_center # radial distance from plate center to injector, to give average distance between injectors
        
        self.angular_distance_injectors = 360/self.N_elements
        self.d_1 = 2*distance_to_center*np.sin(360/self.N_elements/2 * np.pi/180) #right angle triangle between 2 injectors and center of plate, get distance between consecutive injectors

        self.mass_flow_inj = mass_flow_total/N_elements #! further divide by 2 for like doublets

        #* orifice entrance geometry and prefered entrance loss
        self.entrance_type = entrance_type
        self.beta_entrance = beta_entrance
        self.tau_entrance = tau_entrance
        
        #* initial point of study
        self.d_o = d_o
        self.l_o = l_o
        self.l_in = l_in

        #* known fluid properties
        self.gas_density_plenum = gas_density_plenum
        self.gas_R = gas_R
        self.gas_gamma = gas_gamma
        self.gas_T_plenum = gas_T_plenum
        self.gas_dynamic_viscosity = gas_dynamic_viscosity

        # allowed error on convergence of sizing
        self.tolerance = tolerance

        # relative error on d_o, U_jet, C_d
        self.dAlpha = [1, 1]

    def solve (self):

        #* 1) assume value for C_d
        self.C_d = 0.7
        #* 2) calculate dP_inj + gas_pressure_plenum and determine velocity coefficient lambda_2
        self.p_plenum = self.dP_inj + self.p_chamber
        self.lambda_2 = lambda_2(gamma=self.gas_gamma, p_chamber=self.p_chamber, p_plenum=self.p_plenum)
        #* 3) calculate q_lambda_2
        self.q = q_lambda2(gamma=self.gas_gamma, lambda_2=self.lambda_2)
        #* 4) calculate Cstar
        self.Cstar = Cstar(gamma=self.gas_gamma, R=self.gas_R, T_plenum=self.gas_T_plenum)

        # array to record relative error
        self.do_array = [self.d_o]
        self.Cd_array = [self.C_d]
        #
        self.i = 0

        #! iterative sizing
        while max(self.dAlpha)>self.tolerance:
        #* 5) update orifice diameter
            self.d_o = np.sqrt(4/np.pi * self.mass_flow_inj*self.Cstar/(self.C_d*self.p_plenum*self.q))
            self.l_in__d_o = self.l_in/self.d_o

        #* 6) update loss
            # only rounded edge or chamfer possible for gas orifice
            self.xi_in = choose_entrance_loss(entrance_type=self.entrance_type, beta=0, tau__d_o=self.tau_entrance/self.d_o, alpha=0, l_in__d_o=self.l_in__d_o)
            self.d_o__d_1 = self.d_o/self.d_1
            #* 6) compute inlet velocity and Reynolds to estimate friction loss
            self.U_1 = 4/np.pi*self.mass_flow_inj/(self.gas_density_plenum*self.d_o**2) #! ideal injector velocity (at inlet), no contraction
            self.Re = self.gas_density_plenum*self.U_1*self.d_o/self.gas_dynamic_viscosity
            self.friction_coef = blasius_lambda(self.Re)
            self.xi_fric = self.friction_coef*self.l_o/self.d_o

            self.xi = self.xi_in*(1 - self.d_o__d_1**2) + self.xi_fric
            
        #* 7) calculate C_d 
            self.C_d = 1/np.sqrt(1 + self.xi)


            # convergence checks
            self.i+=1
            self.do_array.append(self.d_o)
            self.Cd_array.append(self.C_d)

            # relative error of orifice diameter, Cd and jet velocity
            self.dAlpha = [ abs(self.do_array[self.i]-self.do_array[self.i-1])/self.do_array[self.i-1],\
                            abs(self.Cd_array[self.i]-self.Cd_array[self.i-1])/self.Cd_array[self.i-1]]
        

        self.density_2 = self.gas_density_plenum*(self.p_chamber/self.p_plenum)**(1/self.gas_gamma)
        self.U_2 = gas_exit_velocity(gamma=self.gas_gamma, R=self.gas_R, T_plenum=self.gas_T_plenum, p_chamber=self.p_chamber, p_plenum=self.p_plenum)
        self.momentum = self.density_2*self.U_2**2

#%%





if __name__ == "__main__":

    #* Fuel data 
    # https://iopscience.iop.org/article/10.1088/1742-6596/1473/1/012042/pdf

    p_chamber = 1e5 # max 
    dP__p0 = 0.2
    dP_inj= 1e5 #p_chamber*dP__p0
    mass_flow_tot = 100e-3 # kg/s, approximate total, with +/- 20% margin
    OF_tot = 3.682
    fuel = "gasoline-E10"
    density_fuel = 730 #kg/m3, avg at 20°C, for SP95
    kinematic_viscosity_fuel = 0.5383*1e-6
    oxidizer = "O2"
    massTorch__masstot = 0.0
    mass_flow_torch = massTorch__masstot*mass_flow_tot
    OF_torch = 0.3
    Dh_torch = 50e-3
    Dt_torch = 10e-3
    At_torch = np.pi*Dt_torch**2/4


    mass_flow_main = mass_flow_tot - mass_flow_torch

    mass_flow_ox_main = mass_flow_main*OF_tot/(1 + OF_tot)
    mass_flow_fu_main = mass_flow_main/(1 + OF_tot)

    N_elements = 3
    jets__element = 2

    rhogUg__rhoLUL = 1.5

    d_o=0.5e-3 #! initial orifice diameter
    l_o__d_o=2
    l_in=3e-3
    beta_entrance = 40
    tau_entrance=1e-3
    alpha_entrance=60
    
    # entrance_type="tilted" # "tilted", "chamfer", "round"
    # tolerance=1e-5
    # #* INJECTOR MODEL
    # orifice_flow_test = straight_orifice_flow(  dP_inj=dP_inj, mass_flow_total=mass_flow_fu_main, 
    #                                             N_elements=N_elements, jets__element=jets__element, alpha_entrance=alpha_entrance, entrance_type=entrance_type)
    # orifice_flow_test.solve()
    # print(f"C_d={orifice_flow_test.C_d} | xi={orifice_flow_test.xi}, xi_vortex={orifice_flow_test.xi_1_c}, xi_int={orifice_flow_test.x_in}, xi_fric={orifice_flow_test.xi_fric}")
    # print(f"U_1={orifice_flow_test.U_1} | U_2={orifice_flow_test.U_2} m/s")
    # print(f"dP={dP_inj/1e5} bar")
    # print(f"do={orifice_flow_test.d_o*1e3} | lo={orifice_flow_test.l_o*1e3} mm")
    # print(f"Reynolds entrance={orifice_flow_test.Re}")
    # print(" ")
    # print("gas-liquid jets")
    # print(f"rho_L*U_2^2={orifice_flow_test.momentum} kg/m-s²")
    # print(f"rho_g*Ug^2={orifice_flow_test.momentum*rhogUg__rhoLUL} @ rhogUg/rhoLUL={rhogUg__rhoLUL}")
    # # jet momentum 

