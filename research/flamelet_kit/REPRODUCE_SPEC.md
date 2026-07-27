# Reproduce spec

A from-scratch algorithmic specification, precise enough to re-implement
this kit's two modules (`Flamelet` regime 1, `CoolingPFR` regime 2, plus
`SteadyCache`) without reading the Python. Numbered steps + equations; no
code.

---

## Part A: `Flamelet` (regime 1 -- non-premixed flame in mixture-fraction space)

### A1. Inputs

- A Cantera-loadable reaction mechanism (species thermo/transport/kinetics).
- Number of Z-grid nodes `n_z` (default 65).
- Two boundary streams: oxidizer (`T_ox`, mass-fraction vector `Y_ox`) and
  fuel (`T_fuel`, mass-fraction vector `Y_fuel`), plus an operating pressure `p`.
- A per-species diffusion mask `diff_mask` (default: 1.0 for every species,
  meaning unity Lewis number; entries may be set to 0.0 to make a species
  non-diffusing in Z, Le -> infinity).
- CVODE integrator tolerances `rtol`, `atol` for the chemistry sub-step.

### A2. Stoichiometric mixture fraction (Bilger definition)

For each species k with molecular weight W_k and elemental composition
(n_C(k), n_H(k), n_O(k)) atoms of carbon, hydrogen, oxygen:

    beta(Y) = sum_k  Y_k / W_k * (2*n_C(k) + 0.5*n_H(k) - n_O(k))

    Z_st = ( 0 - beta(Y_ox) ) / ( beta(Y_fuel) - beta(Y_ox) )

Clip Z_st to [0.01, 0.99]. If |beta(Y_fuel) - beta(Y_ox)| < 1e-12 (degenerate
streams), fall back to Z_st = 0.14.

### A3. Z-grid construction (two-sided Roberts/tanh stretching, clustered at Z_st)

Split n_z nodes into a left count `n_left = floor(n_z/2)` and right count
`n_right = n_z - n_left`. Stretching parameter beta = 3.0 (fixed).

Left sub-domain, covering [0, Z_st):
- xi_L = n_left equally spaced points in [0, 1), i.e. `linspace(0,1,n_left+1)`
  with the last point (=1) dropped.
- `Z_L(xi) = Z_st * (1 + sinh(beta*(xi - 1)) / sinh(beta))`
  (fine spacing at xi=1 which maps to Z=Z_st; coarse at xi=0 which maps to Z=0)

Right sub-domain, covering [Z_st, 1]:
- xi_R = n_right equally spaced points in [0, 1] inclusive.
- `Z_R(xi) = Z_st + (1 - Z_st) * sinh(beta*xi) / sinh(beta)`
  (fine spacing at xi=0 which maps to Z=Z_st; coarse at xi=1 which maps to Z=1)

Concatenate Z_L then Z_R into the full grid Z (length n_z, strictly
increasing); pin Z[0]=0.0 and Z[-1]=1.0 exactly to remove floating-point
drift. Verify: for n_z=65, Z_st~0.05-0.5, at least 15 of the 65 nodes fall
within |Z - Z_st| < 0.05.

### A4. Initial condition (linear mixing, frozen/no-reaction)

    T(Z_i)    = T_ox + (T_fuel - T_ox) * Z_i
    Y(Z_i, k) = Y_ox[k]*(1 - Z_i) + Y_fuel[k]*Z_i,     then renormalize each
                row so sum_k Y(Z_i,k) = 1.

Allocate one persistent 0-D constant-pressure reactor per Z-node
(temperature `T(Z_i)` clipped to [250, 4500] K, pressure `p`, composition
`Y(Z_i,:)`), each with its own ODE integrator instance using the supplied
rtol/atol.

### A5. Scalar dissipation profile chi(Z)

Given the caller-supplied scalar `chi_st` (the value of chi at Z_st):

    A(Z)  = erfcinv(2*clip(Z, eps, 1-eps))^2         (eps = 1e-12)
    A_st  = erfcinv(2*clip(Z_st, eps, 1-eps))^2
    chi(Z) = chi_st * exp( 2*A_st - 2*A(Z) )

Note: this profile's mathematical maximum is at Z=0.5 (erfcinv(1)=0), NOT at
Z_st -- chi_st is an anchor value (chi(Z_st) == chi_st exactly), not
necessarily the peak. The diffusion coefficient used everywhere below is
`D(Z) = 0.5 * chi(Z)`.

### A6. Per-step algorithm (Strang splitting)

Given the current state (T, Y arrays over all n_z nodes), a step of size
`dt`, current boundary conditions (`p`, `T_ox`, `Y_ox`, `T_fuel`, `Y_fuel`),
and `chi_st`:

1. Update a pressure-derivative EMA: `dp/dt_raw = (p - p_prev)/dt`;
   `dp_dt_ema = alpha*dp/dt_raw + (1-alpha)*dp_dt_ema`, alpha=0.1 fixed;
   store `p_prev = p`. (This feeds a compression-heating source term below;
   it is identically ~0 if p is constant.)
2. Compute `D(Z) = 0.5*chi(Z)` from A5 using the current `chi_st`.
3. **Diffusion half-step**, `dt/2` (see A7).
4. **Chemistry full-step**, `dt` (see A8).
5. **Diffusion half-step**, `dt/2` (see A7), using the same D(Z) from step 2.
6. Increment step counter and elapsed flamelet time by `dt`.

### A7. Diffusion half-step (Crank-Nicolson, Thomas algorithm)

For each species k with `diff_mask[k] >= 1e-10` (skip diffusion entirely for
species with mask ~0): solve, on the non-uniform Z grid, the 1-D diffusion
PDE `d(phi)/dt = D_k(Z) d^2(phi)/dZ^2` over the half-step `dt_half`, Dirichlet
BCs `phi(Z=0) = Y_ox[k]`, `phi(Z=1) = Y_fuel[k]`, via Crank-Nicolson
(implicit trapezoidal in time) discretized on the non-uniform grid as
follows. Let `h_l(i) = Z(i)-Z(i-1)`, `h_r(i) = Z(i+1)-Z(i)` for interior node
i (1..n_z-2). Second-derivative coefficients at node i:

    a2(i) = 2 / ((h_l+h_r) * h_l)     (coefficient of phi(i-1))
    b2(i) = -2 / (h_l * h_r)          (coefficient of phi(i), negative)
    c2(i) = 2 / ((h_l+h_r) * h_r)     (coefficient of phi(i+1))

Crank-Nicolson gives, for each interior node, an equation of the form
`phi_new(i) - dt_half*D(i)*(a2 phi_new(i-1) + b2 phi_new(i) + c2 phi_new(i+1))
= phi_old(i) + dt_half*D(i)*(a2 phi_old(i-1) + b2 phi_old(i) + c2 phi_old(i+1))`
-- i.e. average the explicit and implicit second-derivative operator (in
this implementation D is evaluated once, outside the Thomas solve, so
effectively the SAME operator appears on both sides scaled by dt_half; this
is the half-explicit/half-implicit CN split). Assemble the tridiagonal
system for the (n_z-2) interior unknowns and solve via the Thomas algorithm
(forward elimination + back substitution). Fold the known Dirichlet
boundary values into the right-hand side of the first and last interior
equations. After solving, clip Y to >= 0 and renormalize each node's Y row
to sum to 1.

Then compute, from the current (post-species-diffusion) state via Cantera at
each node's (T, p, Y): mixture cp_mass, per-species partial molar cp divided
by molecular weight (cp_k), and density rho.

Temperature CN step, over the same `dt_half`, with THREE additional terms
folded in relative to the pure-diffusion species case:
1. An explicit source term `dp_dt_heat(Z) = dp_dt_ema / max(rho(Z)*cp(Z), 1)`
   (compression heating), added to the right-hand side scaled by `dt_half`.
2. A FULLY-IMPLICIT convection term representing the cp-gradient
   ("Pitsch-Peters") flamelet correction:

       dcp/dZ(Z)          via a numerical gradient of cp(Z) over the Z grid
       sum_k cp_k * dY_k/dZ  via a numerical gradient of each species' Y(Z)
       conv(Z) = (D(Z) / max(cp(Z),1)) * (dcp/dZ(Z) + sum_k cp_k(Z)*dY_k/dZ(Z))

   `conv(Z)` is a first-derivative ("advection") coefficient added to the
   SAME tridiagonal system as an additional fully-implicit term (i.e. its
   coefficient is evaluated once at the start of the half-step and held
   fixed through the implicit solve -- this is what "fully implicit, frozen
   coefficient" means; it is NOT an explicit source, it modifies the
   tridiagonal matrix). First-derivative coefficients on the non-uniform
   grid at interior node i:

       a1(i) = -h_r / (h_l*(h_l+h_r))
       b1(i) = (h_r - h_l) / (h_l*h_r)
       c1(i) =  h_l / (h_r*(h_l+h_r))

   fold `-dt_half*conv(i)*a1(i)` into the sub-diagonal, `-dt_half*conv(i)*b1(i)`
   into the diagonal, `-dt_half*conv(i)*c1(i)` into the super-diagonal (same
   sign convention as the diffusion terms), with the boundary-value
   contributions moved to the right-hand side exactly as for diffusion.

   IMPORTANT: this term MUST be fully implicit. An explicit (or
   semi-explicit) treatment was found (in the source rocket solver) to
   overshoot on the Z_st-clustered grid, saturate the temperature clip, and
   drive a spurious high-frequency quench/spike limit cycle. Any
   reimplementation that treats this term explicitly should expect the same
   instability.
3. Dirichlet BCs `T(Z=0) = T_ox`, `T(Z=1) = T_fuel`.

Solve the resulting tridiagonal system via Thomas algorithm, clip result to
[200, 4500] K.

### A8. Chemistry full-step

For each INTERIOR node i (boundary nodes 0 and n_z-1 are Dirichlet and are
never advanced by chemistry): set that node's persistent reactor's state to
(T(i) clipped to [250,4500], p, Y(i,:)), sync the reactor, then integrate the
reactor's ODE system forward by the full `dt` (starting the integrator's
internal clock at 0 each call -- i.e. each node's reactor is advanced by
exactly `dt` of physical time per outer step, independent of prior calls'
internal integrator state, though the THERMOCHEMICAL state itself does
persist node-to-node across steps via the reactor object). Read back the
resulting (T, Y), clip T to [250, 4500], clip Y to >= 0, and store. If the
integrator fails for a node (exception), leave that node's T/Y unchanged for
this step (do not propagate the failure). After the loop, re-impose the
Dirichlet boundary conditions exactly: `T(0)=T_ox`, `T(n_z-1)=T_fuel`,
`Y(0,:)=Y_ox`, `Y(n_z-1,:)=Y_fuel`.

### A9. Outputs / accessors

- `T_max` = max over all nodes of T.
- `T_at_Z(Z_query)`, `Y_at_Z(Z_query, species)`: linear interpolation of the
  T(Z) / Y(Z) arrays onto arbitrary query mixture-fraction value(s) -- this
  is how a host flow field maps its own local Z(x) back to physical T/Y.
- `n_nodes_near_Z_st(half_width)`: count of grid nodes within
  `half_width` of Z_st (grid-quality diagnostic).

---

## Part B: `CoolingPFR` (regime 2 -- single-stream plug-flow reactor with wall heat loss)

### B1. Inputs

- Same Cantera mechanism as regime 1 (for species-index consistency across
  the hand-off).
- Inlet state: `T_in`, mass-fraction vector `Y_in`, pressure `p` (assumed
  constant along the channel), mass flow rate `mdot` (assumed constant along
  the channel -- steady mass conservation, no leakage).
- Channel geometry: `diameter` (hydraulic diameter for a circular
  cross-section, or the appropriate hydraulic diameter otherwise), `length`,
  number of axial segments `n_steps` (so `dx = length/n_steps`).
- Wall model: convective coefficient `h_conv` [W/m^2/K] and wall temperature
  `T_wall` (a constant, or a function of axial position x).

### B2. Per-segment algorithm

Cross-sectional area `A = pi*diameter^2/4`; wetted perimeter `P = pi*diameter`.

Initialize a single persistent 0-D constant-pressure reactor at
`(T_in, p, Y_in)`. For segment i = 0 .. n_steps-1, with running axial
position x (starting at 0):

1. Evaluate `T_wall(x)` (constant or callable).
2. Compute local velocity `u = mdot / (rho(x) * A)` from the reactor's
   CURRENT density.
3. Convert the spatial segment to a time increment: `dt = dx / u`.
4. Chemistry sub-step: advance the (adiabatic) reactor's ODE system forward
   by `dt` (standard constant-pressure ideal-gas reactor integration --
   conserves specific enthalpy exactly in the absence of external heat
   input, since it is adiabatic).
5. Wall heat-loss sub-step: let `T_chem` = the reactor's temperature after
   step 4. Compute the heat-loss rate per unit length
   `q_wall = h_conv * P * (T_chem - T_wall(x))` [W/m], and the segment's
   total heat removed `Q_seg = q_wall * dx` [W]. Read the reactor's current
   specific enthalpy `h_before`. Compute `h_after = h_before - Q_seg/mdot`.
   Set the reactor's state directly via its (enthalpy, pressure, mass
   fraction) state setter to `(h_after, p, Y_now)` -- i.e. re-solve for the
   temperature that gives this new (lower) enthalpy at fixed p and fixed
   composition. This is exact by construction: it does not depend on cp or
   any linearization.
6. Accumulate `Q_wall_total += Q_seg`; advance `x += dx`; record
   `(x, T, Y)` at this segment boundary.

### B3. Energy balance invariant (use to test correctness)

Let `h_in` = the reactor's specific enthalpy at the very start (before
segment 0's chemistry sub-step), and `h_out` = the specific enthalpy after
the final segment's heat-loss sub-step. By construction of B2 step 5,
applied cumulatively:

    mdot * (h_in - h_out)  ==  Q_wall_total     (to reactor/ODE tolerance)

This must hold to a tight relative tolerance regardless of h_conv. As a
special case, with `h_conv = 0` (no wall heat loss), `Q_wall_total = 0` and
`h_out == h_in` (adiabatic reactor conserves enthalpy) -- i.e. temperature at
the outlet approximately equals the inlet temperature (equal exactly if the
reactor is already at its own chemical equilibrium; otherwise it may still
shift somewhat due to ongoing finite-rate reaction, but total enthalpy is
conserved either way).

### B4. Outputs

- Arrays `x`, `T`, `Y` at each of the `n_steps+1` segment boundaries
  (including the inlet).
- `Q_wall_total_W`, `h_in`, `h_out`, `energy_balance_residual_W` (should be
  ~0; see B3).

---

## Part C: `SteadyCache` (cost lever, shared by both regimes)

### C1. State

- Three tolerances: `tol_dT` [same units as the dT you pass in, typically
  K], `tol_p` [relative], `tol_chi` [relative], `tol_T_fuel` [absolute, same
  units as the T_fuel-slot value you pass in].
- Optional EMA time constants `p_ema_tau`, `chi_ema_tau` for filtering the
  cache KEYS only (0 = unfiltered/instantaneous).
- "Last recorded advance" state: `{p, chi, T_fuel, dT}` or "none yet".

### C2. Key computation (optional low-pass filter)

For a raw value `v` (p or chi_st/mdot) and elapsed time `dt_elapsed` since
the last key update: if the corresponding EMA tau is 0, the key is just `v`
(unfiltered). Otherwise maintain an EMA:
`a = 1 - exp(-dt_elapsed/tau)`; `ema = (1-a)*ema_prev + a*v` (initialize
`ema = v` on the first call, or whenever `dt_elapsed <= 0`); the key is the
current `ema`.

### C3. Gate check

Given candidate keys `(p_key, chi_key, T_fuel)`:

1. If there is no prior recorded advance -> MISS, reason "no_prior_advance".
2. Else if the last recorded advance's own `dT` (how much its target
   quantity, e.g. T_max, changed during that advance) exceeds `tol_dT` ->
   MISS, reason "dT". (Rationale: a large last-advance dT means the solution
   was still relaxing, not yet steady -- caching would freeze a transient.)
3. Else if `|p_key - last.p| > tol_p * max(last.p, 1)` -> MISS, reason "p".
4. Else if `|chi_key - last.chi| > tol_chi * max(|last.chi|, 1e-12)` ->
   MISS, reason "chi".
5. Else if `|T_fuel - last.T_fuel| > tol_T_fuel` -> MISS, reason "T_fuel".
6. Else -> HIT (empty-string reason): safe to reuse the existing solution,
   skip re-solving.

### C4. Recording

After actually performing a real advance/march (resolving a MISS), record
the NEW `(p_key, chi_key, T_fuel, dT)` as the "last recorded advance" state
for future gate checks. `dT` here is the change in the solution's own
target diagnostic (e.g. |T_max_after - T_max_before| for a flamelet, or
|T_outlet_after - T_outlet_before| for a PFR) caused by that specific
advance/march.

### C5. Deliberate omission vs. the source

The source rocket cache additionally gated on a "burning" precondition
(upper physical temperature bound, an is-burning check, and a
respark-margin check) before considering any of the above tolerance gates.
That gate exists solely to protect against caching a state on the wrong
branch of an ignition/extinction S-curve. A steady problem (no
ignition/extinction) has no second branch, so this precondition is
correctly omitted here -- do not reintroduce it unless you reintroduce an
ignition/extinction model alongside it.
