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



class straight_orifice_flow_LIQUID ():
    def __init__(self, 
                 dP_inj, p_chamber, mass_flow_total, N_elements, jets__element,
                 #SP95
                 liquid_density=740, liquid_kin_viscosity=0.5383*1e-6,
                 # 
                 d_o=1e-3, l_o=3e-3 , l_in=3e-3,
                 beta_entrance = 40, tau_entrance=1e-3, alpha_entrance=60,
                 entrance_type="tilted", # "tilted", "chamfer", "round"
                 tolerance=1e-6):
        """
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
        self.jets__element = jets__element # if 2 then like doublet, if 3 then triplet

        self.mass_flow_inj = mass_flow_total/N_elements/self.jets__element #! further divide by 2 for like doublets

        #* orifice entrance geometry and prefered entrance loss
        self.entrance_type = entrance_type
        self.beta_entrance = beta_entrance
        self.tau_entrance = tau_entrance
        self.alpha_entrance = alpha_entrance

        #* initial point of study
        self.d_o = d_o
        self.l_o = l_o
        self.l_in = l_in
        self.l_in__d_o = l_in/d_o



        #* known fluid properties
        self.liquid_density = liquid_density
        self.liquid_kin_viscosity = liquid_kin_viscosity
        self.liquid_dyn_viscosity = liquid_kin_viscosity*liquid_density

        # allowed error on convergence of sizing
        self.tolerance = tolerance

        # relative error on d_o, U_jet, C_d
        self.dAlpha = [1, 1, 1]

    def solve (self):

        #* 1) initial data includes :  l_o
        #* 2) begin computation with assumed loss dependent only on entrant loss, xi_in
        

        self.x_in = choose_entrance_loss(entrance_type=self.entrance_type, beta=self.beta_entrance, tau__d_o=self.tau_entrance/self.d_o, alpha=self.alpha_entrance, l_in__d_o=self.l_in__d_o)
        self.C_d = 1/np.sqrt(1 + self.x_in)
        self.U_inj = 4/np.pi*self.mass_flow_inj/(self.liquid_density*self.d_o**2)
        
        # array to record relative error
        self.do_array = [self.d_o]
        self.Cd_array = [self.C_d]
        self.U_array = [self.U_inj]

        #
        self.i = 0

        #! iterative sizing
        while max(self.dAlpha)>self.tolerance:
        #* 3) calculate injector passage diameter, from mass_flow_inj, C_d, liquid_density, dP_inj
            self.d_o = np.sqrt(self.mass_flow_inj*4/np.pi)*np.power(self.C_d, -0.5)*np.power(self.liquid_density*self.dP_inj, -0.25) #discharge coef equation
            self.l_in__d_o = self.l_in/self.d_o

        #* 4) calculate Re and U_inj
            self.U_1 = 4/np.pi*self.mass_flow_inj/(self.liquid_density*self.d_o**2) #! ideal injector velocity (at inlet), no contraction
            self.Re = self.U_1*self.d_o/self.liquid_kin_viscosity
        #* 5) calculate friction loss
            # pipe friction loss
            self.friction_coef = blasius_lambda(self.Re)
            self.xi_fric = self.friction_coef*self.l_o/self.d_o
        #* 6) calculate total loss
            # loss due to inlet vortex
            self.xi_1_c = xi1c_from_fig22(Re=self.Re)
            # entrance loss
            self.x_in = choose_entrance_loss(entrance_type=self.entrance_type, beta=self.beta_entrance, tau__d_o=self.tau_entrance/self.d_o, alpha=self.alpha_entrance, l_in__d_o=self.l_in__d_o)
            # total loss
            self.xi = self.xi_1_c + self.x_in + self.xi_fric
        #* 7) calculate C_d 
            self.C_d = 1/np.sqrt(1 + self.xi)
        #* 8) calculate exit jet velocity, U_2 
            # exit jet velocity, U_2, from dP_inj = density*U_2**2/2(1 + xi)
            self.U_2 = self.C_d * np.sqrt(2*self.dP_inj/self.liquid_density)

        

            # convergence checks
            self.i+=1
            self.do_array.append(self.d_o)
            self.Cd_array.append(self.C_d)
            self.U_array.append(self.U_1)

            # relative error of orifice diameter, Cd and jet velocity
            self.dAlpha = [ abs(self.do_array[self.i]-self.do_array[self.i-1])/self.do_array[self.i-1],\
                            abs(self.Cd_array[self.i]-self.Cd_array[self.i-1])/self.Cd_array[self.i-1],\
                            abs(self.U_array[self.i]-self.U_array[self.i-1])/self.U_array[self.i-1]]
            
        self.eta = (1 + self.xi_1_c + self.x_in + self.xi_fric)/(2*np.sqrt(self.x_in)) #cavitation coefficient
        self.p_n = self.p_chamber - self.dP_inj/self.eta

        self.momentum = self.liquid_density*self.U_2**2

    def solve_mass_flow(self, do):
        """ 
        All conditions fixed (chamber pressure, dP, fluid properties), what is the impact of changing do on the 
        mass flow and C_d ?
        """

        self.do_ = do 
        self.l_in__d_o = self.l_in/self.do_
        
        self.mass_flow_inj_ = self.mass_flow_inj

        self.x_in = choose_entrance_loss(entrance_type=self.entrance_type, beta=self.beta_entrance, tau__d_o=self.tau_entrance/self.do_, alpha=self.alpha_entrance, l_in__d_o=self.l_in__d_o)
        self.C_d_ = self.C_d
        self.U_inj_ = 4/np.pi*self.mass_flow_inj_/(self.liquid_density*self.do_**2)
        
        # array to record relative error
        self.mass_array = [self.mass_flow_inj_]
        self.Cd_array = [self.C_d_]
        self.U_array = [self.U_inj_]

        #
        self.dAlpha = [1, 1, 1]
        self.i = 0

        #! iterative sizing
        while max(self.dAlpha)>self.tolerance:
        #* 3) calculate injector passage diameter, from mass_flow_inj, C_d, liquid_density, dP_inj
            self.mass_flow_inj_ = self.C_d_ * np.pi*self.do_**2 /4 * np.sqrt(self.liquid_density*self.dP_inj)

        #* 4) calculate Re and U_inj
            self.U_1_ = 4/np.pi*self.mass_flow_inj_/(self.liquid_density*self.do_**2) #! ideal injector velocity (at inlet), no contraction
            self.Re_ = self.U_1_*self.do_/self.liquid_kin_viscosity
        #* 5) calculate friction loss
            # pipe friction loss
            self.friction_coef = blasius_lambda(self.Re_)
            self.xi_fric_ = self.friction_coef*self.l_o/self.do_
        #* 6) calculate total loss
            # loss due to inlet vortex
            self.xi_1_c_ = xi1c_from_fig22(Re=self.Re)
            # entrance loss
            self.x_in_ = choose_entrance_loss(entrance_type=self.entrance_type, beta=self.beta_entrance, tau__d_o=self.tau_entrance/self.do_, alpha=self.alpha_entrance, l_in__d_o=self.l_in__d_o)
            # total loss
            self.xi_ = self.xi_1_c_ + self.x_in_ + self.xi_fric_
        #* 7) calculate C_d 
            self.C_d_ = 1/np.sqrt(1 + self.xi_)
        #* 8) calculate exit jet velocity, U_2 
            # exit jet velocity, U_2, from dP_inj = density*U_2**2/2(1 + xi)
            self.U_2_ = self.C_d_ * np.sqrt(2*self.dP_inj/self.liquid_density)

        

            # convergence checks
            self.i+=1
            self.mass_array.append(self.mass_flow_inj_)
            self.Cd_array.append(self.C_d)
            self.U_array.append(self.U_2_)

            # relative error of orifice diameter, Cd and jet velocity
            self.dAlpha = [ abs(self.mass_array[self.i]-self.mass_array[self.i-1])/self.mass_array[self.i-1],\
                            abs(self.Cd_array[self.i]-self.Cd_array[self.i-1])/self.Cd_array[self.i-1],\
                            abs(self.U_array[self.i]-self.U_array[self.i-1])/self.U_array[self.i-1]]


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
    
    entrance_type="tilted" # "tilted", "chamfer", "round"
    tolerance=1e-5
    #* INJECTOR MODEL
    orifice_flow_test = straight_orifice_flow_LIQUID(  dP_inj=dP_inj, mass_flow_total=mass_flow_fu_main, 
                                                N_elements=N_elements, jets__element=jets__element, alpha_entrance=alpha_entrance, entrance_type=entrance_type)
    orifice_flow_test.solve()
    print(f"C_d={orifice_flow_test.C_d} | xi={orifice_flow_test.xi}, xi_vortex={orifice_flow_test.xi_1_c}, xi_int={orifice_flow_test.x_in}, xi_fric={orifice_flow_test.xi_fric}")
    print(f"U_1={orifice_flow_test.U_1} | U_2={orifice_flow_test.U_2} m/s")
    print(f"dP={dP_inj/1e5} bar")
    print(f"do={orifice_flow_test.d_o*1e3} | lo={orifice_flow_test.l_o*1e3} mm")
    print(f"Reynolds entrance={orifice_flow_test.Re}")
    print(" ")
    print("gas-liquid jets")
    print(f"rho_L*U_2^2={orifice_flow_test.momentum} kg/m-s²")
    print(f"rho_g*Ug^2={orifice_flow_test.momentum*rhogUg__rhoLUL} @ rhogUg/rhoLUL={rhogUg__rhoLUL}")
    # jet momentum 

