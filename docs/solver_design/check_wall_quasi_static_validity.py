"""
Wall-model guardrail script backing DESIGN_PLAN_shellntube_transient.md section 4.1.
Also serves as WP3 acceptance test #2 (see the design doc's WP3 validation list).

QUESTION
--------
For the transient solver, how should each axial node's wall be modelled during a
fast He mass-flow ramp (0 -> full in t_ramp seconds, which drives h_c up over ~3
orders of magnitude)? Three candidate models are compared against a
fully-resolved reference:

  A  (truth)   Resolved radial 1D transient conduction across the wall thickness,
               implicit finite-difference, Robin BCs both faces. Grid/timestep
               converged -- the ground truth.

  B  (chosen)  LUMPED wall-energy ODE:  rho*cp*delta * dTbar/dt = q_hot - q_cold,
               with the two face temperatures reconstructed each step from a
               QUADRATIC quasi-static profile consistent with both Robin BCs and
               the current mean. One ODE state per axial node. This is the
               architecture the design doc adopts.

  C  (strawman) Zero-thermal-mass steady resistance network: wall temperature
               teleports to its new steady value the instant h_c changes. No
               inertia. This was mistakenly tested as if it were "the quasi-static
               model" in an earlier draft and wrongly used to reject model B.

RESULT (t_ramp=2.5 s; INCO718, delta=0.85 mm, h_g=300, h_c 300->75000 W/m2K)
---------------------------------------------------------------------------
  A: reference. Wall dT peaks ~68.9 K at t~0.34 s (5x the 13.8 K steady value)
     -- this transient overshoot is the thermal-shock signal we must capture.
  B: max face-T error vs A ~ 1.8 K, flux error ~1-2 %. Captures the 69 K peak
     essentially exactly. --> VALIDATED, this is the architecture.
  C: max face-T error vs A ~ 438 K. Misses the physics entirely. --> rejected.

So the lumped-ODE + quadratic-reconstruction model (B) is both cheap (1 state
per node, a 2x2 solve per step) and accurate (<2 K) for THIN walls. No radial
sub-grid is needed. The design doc's section 4.1 fallback (switch to a radial
sub-grid) only triggers if B's error exceeds threshold for a future config with
thick thermal mass (e.g. an 8 mm shell wall: delta^2/(pi^2 alpha) ~ 11 s, profile
no longer quasi-static) or sub-100 ms ramps.

IMPLEMENTATION WARNING: the 2x2 face-reconstruction algebra (model B) is easy to
get subtly wrong -- a sign error in the cold-face constraint produced diverging
garbage (wall dT > 1000 K) before it was caught. The steady-state
self-consistency assertion below (reconstructed faces == steady-network faces at
fixed h) is the guard; keep an equivalent assertion in the production code.

Re-run with the base interpreter (the repo .venv is broken -- stale path):
  C:\\Users\\raphael.aubry\\AppData\\Local\\Programs\\Python\\Python313\\python.exe
Tweak h_g / h_c_full / t_ramp / material to re-check the B-vs-A gate for any new
config or ramp before trusting the lumped model there.
"""
import numpy as np

# --- material: Inconel 718, representative mid-range values ---
k_w   = 20.0      # W/m-K
rho_w = 8200.0    # kg/m3
cp_w  = 435.0     # J/kg-K   (RT value; c_p(T) table to be added in WP3a, section 4.7)
alpha = k_w/(rho_w*cp_w)
delta = 0.85e-3   # m, wall thickness

# --- boundary conditions ---
T_g = 1500.0      # K, representative bulk hot gas (fixed -- isolating the He-ramp effect)
h_g = 300.0       # W/m2K, fixed (hot side not ramping in this test)
T_c = 400.0       # K, representative bulk He (fixed)
h_c_low  = 300.0     # W/m2K, floor at ~zero/natural-convection He flow
h_c_full = 75000.0   # W/m2K, full-flow value (from a converged run of main_solve.py)
t_ramp = 2.5         # s -- the intended He mass-flow ramp duration

def h_c_of_t(t):
    if t <= 0: return h_c_low
    if t >= t_ramp: return h_c_full
    return h_c_low + (h_c_full - h_c_low)*(t/t_ramp)

# --- steady 1D series-resistance network (used by model C, and as B's self-test target) ---
def steady_solve(h_c):
    R_g = 1.0/h_g
    R_w = delta/k_w
    R_c = 1.0/h_c
    q = (T_g - T_c)/(R_g + R_w + R_c)    # W/m2
    T_wg = T_g - q*R_g
    T_wc = T_wg - q*R_w
    return q, T_wg, T_wc

# =========================================================================
# Model A: resolved radial PDE (truth), backward Euler, half-cell Robin BCs
# =========================================================================
N = 41
dx = delta/(N-1)
dt = 1e-4
q0, Twg0, Twc0 = steady_solve(h_c_low)
x = np.linspace(0, delta, N)
T_field = Twg0 - (Twg0-Twc0)*(x/delta)   # linear steady IC
r = alpha*dt/dx**2
half = rho_w*cp_w*(dx/2)/dt

def stepA(T, h_c):
    Nn = len(T)
    A = np.zeros((Nn, Nn)); b = np.zeros(Nn)
    for i in range(1, Nn-1):
        A[i, i-1] = -r; A[i, i] = 1+2*r; A[i, i+1] = -r; b[i] = T[i]
    A[0, 0] = half + h_g + k_w/dx; A[0, 1] = -k_w/dx; b[0] = half*T[0] + h_g*T_g
    A[-1, -1] = half + h_c + k_w/dx; A[-1, -2] = -k_w/dx; b[-1] = half*T[-1] + h_c*T_c
    return np.linalg.solve(A, b)

# =========================================================================
# Model B: lumped energy ODE + quadratic quasi-static face reconstruction
# =========================================================================
# Profile: T(x) = T_wg - (q_h/k)x + (q_h - q_c) x^2 / (2 k delta),
#   q_h = h_g (T_g - T_wg),   q_c = h_c (T_wc - T_c)
# Constraint 1 (mean == Tbar):  Tbar = T_wg - q_h*(delta/2k) + (q_h-q_c)*(delta/6k)
# Constraint 2 (cold-face value identity from the profile at x=delta):
#   T_wc = T_wg - q_h*(delta/k) + (q_h-q_c)*(delta/2k)
#        = T_wg - (q_h+q_c)*delta/(2k)
# Linear 2x2 in (T_wg, T_wc):
Aq = delta/(2*k_w); Bq = delta/(6*k_w); Cq = delta/k_w
def faces_from_Tbar(Tbar, h_c):
    M = np.array([[1 + h_g*(Aq - Bq),   -h_c*Bq        ],
                  [-(1 + h_g*Cq/2),      1 + h_c*Cq/2  ]])
    rhs = np.array([Tbar + h_g*T_g*(Aq - Bq) - h_c*T_c*Bq,
                    h_c*T_c*Cq/2 - h_g*T_g*Cq/2])
    T_wg, T_wc = np.linalg.solve(M, rhs)
    q_h = h_g*(T_g - T_wg); q_c = h_c*(T_wc - T_c)
    return T_wg, T_wc, q_h, q_c

# ---- self-consistency assertion: at fixed h, reconstruction == steady network ----
qs, TwgS, TwcS = steady_solve(h_c_full)
TbarS = TwgS - qs*delta/(2*k_w)   # mean of the linear steady profile
w_chk, c_chk, _, _ = faces_from_Tbar(TbarS, h_c_full)
assert abs(w_chk - TwgS) < 1e-6 and abs(c_chk - TwcS) < 1e-6, \
    f"B face-reconstruction FAILS steady self-test: {w_chk:.4f}/{c_chk:.4f} vs {TwgS:.4f}/{TwcS:.4f}"

Tbar = Twg0 - q0*delta/(2*k_w)    # lumped IC consistent with steady low-flow profile

# =========================================================================
# March all three models on the same clock
# =========================================================================
n_steps = int(round(t_ramp*1.6/dt))
t = 0.0
errB_wg = errB_wc = errB_q = 0.0; tB = 0.0
errC_wg = errC_wc = 0.0
maxdT_A = 0.0; tdT = 0.0; dTB_at_peak = 0.0
samples = [0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0]; si = 0

print("   t      h_c    | A(truth) Twg   Twc   dT | B(lumped) Twg   Twc   dT | C(strawman) Twg   Twc")
for _ in range(n_steps+1):
    h_c = h_c_of_t(t)
    TwgA, TwcA = T_field[0], T_field[-1]; dTA = TwgA - TwcA
    TwgB, TwcB, qhB, qcB = faces_from_Tbar(Tbar, h_c); dTB = TwgB - TwcB
    _, TwgC, TwcC = steady_solve(h_c)
    if t <= t_ramp*1.3:
        if abs(TwgA-TwgB) > errB_wg: errB_wg = abs(TwgA-TwgB); tB = t
        errB_wc = max(errB_wc, abs(TwcA-TwcB))
        errB_q = max(errB_q, abs(h_c*(TwcA-T_c) - qcB))
        errC_wg = max(errC_wg, abs(TwgA-TwgC)); errC_wc = max(errC_wc, abs(TwcA-TwcC))
        if dTA > maxdT_A: maxdT_A = dTA; tdT = t; dTB_at_peak = dTB
    if si < len(samples) and abs(t-samples[si]) < dt/2:
        print(f"{t:6.3f} {h_c:8.0f} |   {TwgA:7.1f} {TwcA:7.1f} {dTA:4.1f} |    {TwgB:7.1f} {TwcB:7.1f} {dTB:4.1f} |     {TwgC:7.1f} {TwcC:7.1f}")
        si += 1
    T_field = stepA(T_field, h_c)
    Tbar = Tbar + (qhB - qcB)/(rho_w*cp_w*delta)*dt
    t += dt

print()
print("Steady self-consistency assertion (B faces == steady network): PASS")
print(f"B  (lumped ODE + quad profile)  vs truth: max|Twg|={errB_wg:.2f} K @t={tB:.3f}s, "
      f"max|Twc|={errB_wc:.2f} K, max|q_c|={errB_q/1e3:.1f} kW/m2   -> {'PASS' if max(errB_wg,errB_wc)<5 else 'FAIL'} (<5 K gate)")
print(f"C  (zero-mass steady network)   vs truth: max|Twg|={errC_wg:.1f} K, max|Twc|={errC_wc:.1f} K   -> rejected")
print(f"Transient wall-dT peak (truth) = {maxdT_A:.1f} K @t={tdT:.3f}s "
      f"(B reproduces {dTB_at_peak:.1f} K); steady-state dT = {steady_solve(h_c_full)[1]-steady_solve(h_c_full)[2]:.1f} K")
