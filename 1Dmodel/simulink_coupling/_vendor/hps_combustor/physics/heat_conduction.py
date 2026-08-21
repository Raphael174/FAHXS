""" 
@author : Raphaël Aubry

Functions to compute the conductance and heat conduction resolution in 1D steady across different
layers of material

UP = U*P where U=convection coefficient and P=wetted perimeter
    so UP with dx = 1 gives UA which is the actual correct definition of conductance
    
UP is used for the 1D variation of the heat flux along the length of the engine(channel)

"""

# Library Imports
import numpy as np
from scipy.optimize import fsolve
from scipy.special import i0, iv  # Modified Bessel function of the first kind - for eta_fin of triangular profile
# --- patched snippet inside your class file ---
from .radiation_model.radiation_equations import qrad_net_mbl, hrad_from_q
######## CONDUCTANCE FOR VARIOUS LAYER COMBINATIONS 


def compute_UP_singlewall(h_g, h_c, P_g, P_c, P_w, s_w, k_w, h_g_rad=0, dx=1):
    """ 
    Returns perimeter-wise conductance from hot gas to coolant and through a single wall
    
    Evaluate wall k_w according to average temperature from hot to cold side
    
    P_g, P_w, P_c : Wetted perimeter at given radial distance from center
                    P_g and P_w have the same value since they are taken at the hot gas wall
                    P_c is cummulative perimeter of cooling channels (when multiplied by dx you get the coolant wall wet surface area)
                    
                    
            The value of the P_g and P_w shall be the total value divided by the number of channels
            If not we are basically saying that all heat is focusing on a single channel
            
            
    """
    return (1/((h_g + h_g_rad)*P_g*dx) + s_w/(k_w*P_w*dx) + 1/(h_c*P_c*dx))**(-1)



def compute_fin_efficiency_ch (h_c, k_w, t_w, h_ch):
    """ 
    fin efficiency assuming adiabatic tip and very long fin
    """
    m = np.sqrt(h_c/(k_w*t_w)) #! factor x1 since one fluid touching 1 side of each fin
    eta_fin = np.tanh(m*h_ch)/(m*h_ch)

    return eta_fin, m

def compute_fin_efficiency_Pizzarelli (h_c, k_w, t_w, h_ch):
    """ 
    fin efficiency assuming adiabatic tip and very long fin
    """
    m = np.sqrt(2*h_c/(k_w*t_w)) #! factor x1 since one fluid touching 1 side of each fin
    eta_fin = np.sqrt(2*k_w/(h_c*t_w)) * np.tanh(m*h_ch)

    return eta_fin, m

"""
from scipy.special import iv  # Modified Bessel function of the first kind - for eta_fin of triangular profile
 if fin_profile == "rectangular":
            A_s_fin = 2 * w_fin * L_fin # For adiabatic tip, for convective tip add A_c_tip #!  #! Valid for adiabatic tip and constant cross section
            eta_fin = np.tanh(mL) / mL # Fin efficiency  #! Valid for adiabatic tip and constant cross section
            #q_dot_fin = (T_wh - T_h) * np.sqrt(h_h * per_fin * k_fin * A_fin_base) * np.tanh(mL) #! Valid for adiabatic tip and constant cross section
        else:
            A_s_fin = 2 * w_fin * np.sqrt(L_fin**2 + (t_fin_base/2)**2)
            eta_fin = iv(1, 2*mL) /(mL * iv(0, 2*mL))
 """

def compute_eta_fin_rectangular (h_g, k_fin, t_fin, height_fin):
    """Adiabatic tip, constant rectangular cross section"""
    m = np.sqrt(2*h_g/(k_fin*t_fin)) #! factor x2 since one fluid touching 2 sides of fin
    eta_fin = np.tanh(m*height_fin)/(m*height_fin)
    
    return eta_fin, m

def compute_eta_fin_triangular (h_g, k_fin, t_fin, height_fin):
    """Adiabatic tip, constant rectangular cross section"""
    m = np.sqrt(2*h_g/(k_fin*t_fin))
    eta_fin = iv(1, 2*m*height_fin) /(m*height_fin * i0(2*m*height_fin))
    
    return eta_fin, m


#%%


class OneDimensionalSteadyConduction_ShellnHelicalTube : 
    
    """ 
    1D planar wall conduction
    Written as class to leverage the iterative solve of wall temperatures and thermal conductivity
    """
    
    def __init__(self,
                    h_g, h_c,
                    T_c, T_g, 
                    s_w, Dh_ch,
                    f_kw_at_T, 
                    T_wg_0, T_wc_0, T_c_check_0, 
                    dx=1,
                    # --- new radiation knobs ---
                    rad_enabled=False, eps_s=0.80,
                    rad_backend=None,   # callable: (T_eval, state_dict)-> epsilon_g(T_eval)
                    rad_state=None,     # dict holding {'p':..., 'yH2O':..., 'yCO2':..., 'Le':...}
                    # "outer" (default): hot fluid outside the tube, cold inside — the
                    # helical-coil config (combustion gas in shell, He in coil).
                    # "inner": hot fluid INSIDE the tube, cold outside — the shell-and-tube
                    # config (combustion gas in tubes, coolant on the shell side). Only
                    # affects which perimeter (inner πDh vs outer π(Dh+2s)) each flux uses;
                    # the thin-wall quadratic profile factors (a2,a6) are unchanged since
                    # they depend only on thickness/conductivity, not orientation.
                    hot_side="outer",
                    ):
        self.rad_enabled = rad_enabled
        self.eps_s = eps_s
        self.rad_backend = rad_backend
        self.rad_state = rad_state or {}
        self.hot_side = hot_side


        self.h_g = h_g 
        self.h_c = h_c 

        self.T_g = T_g 
        self.T_c = T_c 
            
        self.s_w = s_w
        self.Dh_ch = Dh_ch

        # function to determine thermal conductivity of wall wrt T
        self.f_kw_at_T = f_kw_at_T 
        
        # initial guess wall temperatures 
        self.T_wg_0 = T_wg_0 
        self.T_wc_0 = T_wc_0 
        self.T_c_check_0 = T_c_check_0

        self.dx = dx
        if self.dx<=0.:
            self.dx = 1e-6

    
    def Solve1Dconduction(self):
                
        def f_x(temp):
            
            T_wg, T_wc, T_wc_check = temp[0], temp[1], temp[2]

            # average temperatures for thermal conductivity of material
            self.T_w_avg = (T_wg+T_wc)/2 
            # average thermal conductivities
            self.k_w = self.f_kw_at_T(self.T_w_avg)
        
            # --- geometry/areas (per your notation) ---
            P_inner = np.pi * self.Dh_ch
            P_outer = np.pi * (self.Dh_ch + 2 * self.s_w)
            if self.hot_side == "inner":
                P_hot, P_cold = P_inner, P_outer
            else:
                P_hot, P_cold = P_outer, P_inner
            A_hot = P_hot * self.dx
            A_cold = P_cold * self.dx
            
            # starting here from the resistance R and integrating to R*Area
            self.R_wall = np.log((self.Dh_ch/2+self.s_w)/(self.Dh_ch/2))/(2*np.pi*self.dx*self.k_w)
            self.Rdx_wall = self.R_wall*self.dx
            self.RPdx_wall = self.Rdx_wall * P_cold


            # --- radiation: build h_rad on the hot side and add in parallel with h_g ---
            if self.rad_enabled and (self.rad_backend is not None):
                # gas emissivity at Tg (emission) and absorptivity at wall temp (≈ emissivity at Ts)
                self.eps_emit = self.rad_backend(T_eval=self.T_g,  **self.rad_state)  # ε_g(Tg)
                self.eps_abs  = self.rad_backend(T_eval=T_wg,      **self.rad_state)  # α_g(Twg)≈ε_g(Twg)
                self.q_w_rad  = qrad_net_mbl(self.T_g, T_wg, self.eps_emit, self.eps_abs, self.eps_s)  # W/m^2
                self.h_g_rad    = hrad_from_q(self.T_g, T_wg, self.q_w_rad)                         # W/m^2-K
            else:
                self.h_g_rad, self.q_w_rad, self.eps_emit, self.eps_abs= 0., 0., 0., 0.

            self.h_g_eff = self.h_g + self.h_g_rad  # parallel paths

            # --- your U/UA forms with h_g replaced by h_g_eff ---
            self.U  = (1/self.h_g_eff + self.RPdx_wall + 1/(self.h_c))**(-1)
            self.UP = (1/(self.h_g_eff*P_hot) + self.Rdx_wall + 1/(self.h_c*P_cold))**(-1)
            self.UA = (1/(self.h_g_eff*A_hot) + self.R_wall + 1/(self.h_c*A_cold))**(-1)

            # --- heat transfer and wall temps (unchanged algebra) ---
            self.dq__dx = (self.T_g - self.T_c) * self.UP
            self.dQ     = (self.T_g - self.T_c) * self.UA
            self.q_w    = self.dQ / A_hot  # inner-side flux

            self.Res_g = 1/(self.h_g_eff*A_hot)
            self.Res_c = 1/(self.h_c*A_cold)
            self.Res_w = np.copy(self.R_wall)

            self.T_wg_new = self.T_g - self.dQ / (self.h_g_eff*A_hot)
            self.T_wc_new = self.T_wg_new - self.dQ * self.R_wall
            self.T_c_check = self.T_wc_new - self.dQ / (self.h_c*A_cold)

            return [T_wg - self.T_wg_new, T_wc - self.T_wc_new, T_wc_check - self.T_c_check]        
        
        # solving for the wall temperatures based on converged thermal conductivity of wall
        self.T_wg, self.T_wc, self.T_c_check_f = fsolve(func=f_x, x0=[self.T_wg_0, self.T_wc_0, self.T_c_check_0], xtol=1e-8)
        # back chaking that the solution has converged by re-computing coolant temperature from coolant wall


        return self.T_wg, self.T_wc

    # ------------------------------------------------------------------
    # TRANSIENT support: flux evaluation at a prescribed mean wall temp.
    # ------------------------------------------------------------------
    def _faces_from_hgeff(self, T_bar, h_g_eff, a2, a6):
        """Single 2x2 solve for (T_wg, T_wc) at a fixed effective hot-side coefficient.
        Closed form (no np.linalg.solve — this is on the transient hot path)."""
        h_c = self.h_c
        # M = [[m00, m01],[m10, m11]], rhs = [r0, r1]
        m00 = 1.0 + h_g_eff * (a2 - a6); m01 = -h_c * a6
        m10 = -(1.0 + h_g_eff * a2);      m11 = 1.0 + h_c * a2
        r0 = T_bar + h_g_eff * self.T_g * (a2 - a6) - h_c * self.T_c * a6
        r1 = h_c * self.T_c * a2 - h_g_eff * self.T_g * a2
        det = m00 * m11 - m01 * m10
        T_wg = (r0 * m11 - m01 * r1) / det
        T_wc = (m00 * r1 - r0 * m10) / det
        return T_wg, T_wc

    def fluxes_at_Tbar(self, T_bar, h_g_rad=None):
        """
        Transient-solver companion to Solve1Dconduction().

        Given the *lumped* thickness-mean wall temperature T_bar (the single
        time-integrated state per axial node, see
        DESIGN_PLAN_shellntube_transient.md section 4.1), reconstruct the two
        face temperatures from a quadratic quasi-static radial profile and
        return the per-unit-length heat fluxes on each face. The transient
        driver integrates  (rho*cp*A_wall) dT_bar/dt = dq_hot__dx - dq_cold__dx.

        Model B in doc/check_wall_quasi_static_validity.py — validated to <2 K
        against a fully-resolved radial PDE across a fast He ramp.

        Parameters
        ----------
        h_g_rad : float or None
            FAST PATH (transient default): if the caller has already computed the
            radiation coefficient at its current T_wg estimate (e.g. from a
            tabulated emissivity, and owning the outer T_wg fixed point itself),
            pass it here — this does a single closed-form face solve with
            h_g_eff = h_g + h_g_rad and no radiation-backend calls.
            SLOW PATH (h_g_rad=None): iterate the radiation backend internally on
            T_wg, as the steady node does — used for standalone/validation calls.

        Returns dict:
            dq_hot__dx, dq_cold__dx : [W/m]  (hot into wall, cold out of wall)
            T_wg, T_wc              : [K] face temperatures
            h_g_rad, q_w_rad        : radiation diagnostics (0 if disabled)
            k_w                     : wall conductivity at the mean temperature

        Steady self-consistency: at a fixed operating point the fluxes vanish to
        equality (dq_hot=dq_cold) and the reconstructed faces reproduce the
        steady 3-resistance network — the guardrail script asserts this.
        """
        s, Dh = self.s_w, self.Dh_ch
        P_inner = np.pi * Dh
        P_outer = np.pi * (Dh + 2.0 * s)
        if self.hot_side == "inner":
            P_h, P_c = P_inner, P_outer    # hot fluid inside the tube (shell-and-tube config)
        else:
            P_h, P_c = P_outer, P_inner    # hot fluid outside the tube (helical-coil config)
        k_w = self.f_kw_at_T(T_bar)
        a2 = s / (2.0 * k_w)           # δ/2k
        a6 = s / (6.0 * k_w)           # δ/6k

        if h_g_rad is not None:
            # --- fast path: radiation coefficient supplied, single closed-form solve ---
            h_g_eff = self.h_g + h_g_rad
            T_wg, T_wc = self._faces_from_hgeff(T_bar, h_g_eff, a2, a6)
            q_w_rad = h_g_rad * (self.T_g - T_wg)
        else:
            # --- slow path: iterate radiation backend on T_wg internally ---
            T_wg = self.T_wg_0
            h_g_rad, q_w_rad = 0.0, 0.0
            for _ in range(4):
                if self.rad_enabled and (self.rad_backend is not None):
                    eps_emit = self.rad_backend(T_eval=self.T_g, **self.rad_state)
                    eps_abs = self.rad_backend(T_eval=T_wg, **self.rad_state)
                    q_w_rad = qrad_net_mbl(self.T_g, T_wg, eps_emit, eps_abs, self.eps_s)
                    h_g_rad = hrad_from_q(self.T_g, T_wg, q_w_rad)
                else:
                    h_g_rad, q_w_rad = 0.0, 0.0
                h_g_eff = self.h_g + h_g_rad
                T_wg_new, T_wc = self._faces_from_hgeff(T_bar, h_g_eff, a2, a6)
                if abs(T_wg_new - T_wg) < 1e-9:
                    T_wg = T_wg_new
                    break
                T_wg = T_wg_new
            h_g_eff = self.h_g + h_g_rad

        dq_hot__dx = h_g_eff * P_h * (self.T_g - T_wg)
        dq_cold__dx = self.h_c * P_c * (T_wc - self.T_c)

        return dict(dq_hot__dx=dq_hot__dx, dq_cold__dx=dq_cold__dx,
                    T_wg=float(T_wg), T_wc=float(T_wc),
                    h_g_rad=float(h_g_rad), q_w_rad=float(q_w_rad), k_w=float(k_w))
        
