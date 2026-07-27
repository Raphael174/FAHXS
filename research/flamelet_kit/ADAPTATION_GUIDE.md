# Adaptation guide

How to retarget `flamelet_kit` to your own two-regime (flame + cooling
channel) problem.

## 1. Mechanism + fuel/oxidizer streams (regime 1: `Flamelet`)

Swap the Cantera mechanism and stream definitions passed to `Flamelet`:

```python
fl = Flamelet("your_mechanism.yaml", n_z=65)
fl.init_mixing(T_ox=..., Y_ox=your_oxidizer_Y, T_fuel=..., Y_fuel=your_fuel_Y, p=...)
```

**Diesel note**: Cantera does not bundle a diesel mechanism. Diesel is
conventionally represented by a surrogate fuel -- commonly n-dodecane
(n-C12H26), sometimes blended with n-decane / n-heptane / a small aromatic
fraction for a closer cetane-number/sooting match. Practical path:
1. Obtain a reduced n-dodecane (or your chosen surrogate) Cantera-format
   mechanism -- e.g. a reduced LLNL n-dodecane mechanism converted to YAML
   (search "n-dodecane reduced mechanism Cantera yaml"; several
   35-100-species reduced mechanisms are published and commonly converted
   with `ck2yaml`/`cantera.ck2yaml` from CHEMKIN format).
2. Point `Flamelet("path/to/your_dodecane_mechanism.yaml")` at it.
3. Define `Y_fuel` as pure (or your blend's) surrogate species mass
   fractions, `Y_ox` as your oxidizer (pure O2, or O2 diluted with an inert
   if that matches your actual injector).
4. Re-derive `Z_st` automatically (the Bilger formula in `flamelet.py` needs
   only the mechanism's species compositions -- nothing diesel-specific).

## 2. Grid resolution (n_z)

Default 65 nodes. More nodes = finer resolution near Z_st at higher
per-step chemistry cost (chemistry integrates one Cantera reactor per
interior node). 31-41 is adequate for quick iteration/prototyping; use 65+
for a case you intend to trust quantitatively, and always confirm you have
converged with a simple n_z doubling check before trusting quench/burn
outcomes at a given chi_st.

## 3. Supplying chi_st (regime 1)

The source rocket solver derived chi_st from a wall-shear turbulent-strain
closure specific to the grain boundary layer -- not portable. For a heat
exchanger's mixing/flame zone, three concrete estimation routes (pick
whichever matches what your outer flow model actually resolves):

1. **Gradient definition** (most direct, if you resolve/track Z as a
   transported scalar): `chi = 2*D*|grad(Z)|^2` at the stoichiometric
   surface, D = mass diffusivity (or turbulent diffusivity if using a
   turbulence model). Evaluate |grad Z| numerically from your CFD/1-D
   solution near Z_st.
2. **Scalar mixing frequency** (if you have a turbulence closure with
   turbulent kinetic energy k and dissipation epsilon): `chi ~ C * epsilon/k`
   for an O(1) constant C (commonly ~2, tune/calibrate against a known
   case) -- the standard turbulent scalar-mixing-frequency estimate.
3. **Residence/strain scaling** (simplest, if you know only bulk flow
   parameters): `chi_st ~ U / L`, U = characteristic injection/mixing
   velocity, L = characteristic mixing-layer thickness (jet width, injector
   diameter, or similar). This is a scaling estimate, not a closure -- use it
   to get the right order of magnitude, then refine with (1) or (2) if you
   need quantitative accuracy.

Whichever route you choose, `chi_st` is just an argument to `Flamelet.step`
-- there is no other coupling point to touch.

## 4. Boundary stream states

`T_ox`/`Y_ox` and `T_fuel`/`Y_fuel` are passed to every `step()` call (they
can change call-to-call if your feed conditions ramp) -- there is no
persistent stream-state object to update elsewhere.

## 5. Cache tolerances (`SteadyCache`, both regimes)

Defaults (`tol_dT=2.0` K, `tol_p=0.01` rel, `tol_chi=0.05` rel,
`tol_T_fuel=5.0` K) came from the rocket source and are a reasonable
starting point. Tune them against YOUR transient rate:
- If your operating point ramps slowly (minutes), tighten tolerances (cache
  is nearly always safe to hit; you want to catch real drift).
- If you have any high-frequency ripple you don't want to chase (e.g. small
  numerical noise in an upstream flow solve), either widen tolerances
  slightly or use the optional EMA-filtered keys (`p_ema_tau`, `chi_ema_tau`
  on `SteadyCache`) to low-pass the cache KEY only -- this changes only the
  skip/recompute decision, never what the flamelet/PFR actually integrates
  with on a miss.
- Watch `cache.hit_rate` / `cache.miss_reasons` in your own instrumentation;
  a low hit rate with most misses on one gate tells you exactly which input
  is moving too fast for the current tolerance.

## 6. Regime-2 (`CoolingPFR`) retargeting

- Reuses the SAME mechanism as regime 1 (species indices must line up for
  the T_at_Z/Y_at_Z hand-off to make sense).
- Channel geometry: `diameter` (hydraulic diameter for non-circular ducts),
  `length`, `n_steps` (axial resolution).
- Heat loss: `h_conv` (a convective coefficient you supply -- from your own
  wall correlation, e.g. Dittus-Boelter/Nusselt correlation for the channel
  Reynolds number) and `T_wall` (constant, or `callable(x) -> T_wall(x)` for
  an axially-varying wall temperature, e.g. from a coupled wall conduction
  solve).
- Caching: `CoolingPFR.march_cached(cache, ...)` re-purposes `SteadyCache`'s
  three generic tolerance slots (p/chi/T_fuel) as (p, mdot, T_in) -- see the
  docstring in `cooling_pfr.py` for the exact mapping.

## 7. WHAT WAS REMOVED AND WHY

Everything below existed in the source rocket RIF solver to handle
ignition and extinction on a transient S-curve. A steady or
steadily-changing heat exchanger with an already-established flame (regime
1) feeding an already-flowing channel (regime 2) has no such transient to
track, so all of it was deliberately stripped:

- **`spark()` / HP-equilibrium re-ignition**: re-initializes nodes to their
  adiabatic-equilibrium state to force-ignite a cold flamelet. Not needed --
  a steady flame is assumed already burning; if you need to seed a new
  representative flamelet from cold, do it once, out-of-band, the way
  `example_run.py::seed_flame_kernel` demonstrates (explicitly NOT part of
  the library, since it's a one-time initialization convenience, not a
  runtime ignition model).
- **Quench/extinction predicate** (`is_burning`, minimum-burning-temperature
  floor, consecutive-quench counting, respark margins): exists to detect
  when the flamelet has dropped off the flame branch onto the quenched
  branch of the S-curve, and to gate re-ignition/caching around that. A
  steady flow has (by the problem's own scoping) no second branch to fall
  onto.
- **Reactive-species hot-band chemistry gating**: an optimization that only
  advances chemistry near an active ignition/flame front to save cost during
  a transient. This kit always advances chemistry at every interior node --
  simpler and correct for a steady problem where the flame position doesn't
  move around.
- **NN-surrogate shadow evaluation, transition recorder, all telemetry
  hooks**: infrastructure for a separate machine-learning-surrogate
  validation campaign specific to the rocket project. Entirely orthogonal to
  the numerics and not needed to reproduce or use the flamelet/PFR solvers.
- **CSOLID/BIN soot species diff_mask special-casing**: `diff_mask` is now a
  clean, generic optional constructor argument defaulting to all-ones
  (unity Lewis number for every species). Supply your own mask if your
  mechanism has a species you want to treat as non-diffusing in Z (Le -> inf).

## 8. Applicability caveat (read before committing to this approach)

The mixture-fraction flamelet abstraction (`Flamelet`, regime 1) is valid
ONLY for a **non-premixed, two-feed-stream, diffusion-controlled** reacting
flow, where a conserved scalar Z is physically meaningful (it must be 0 in
one pure feed stream, 1 in the other, and materially conserved -- no source
term -- everywhere between). A diesel/O2 injector where fuel and oxidizer
arrive separately fits this description and is squarely in the valid domain.

If your regime-1 mixing/reaction zone is instead **premixed** (fuel and
oxidizer already mixed before reacting) or **single-stream** (only one feed,
e.g. a pre-vaporized/pre-mixed charge), the flamelet-in-Z model is the WRONG
reduced model -- there is no second stream to define Z against, and forcing
one in will not represent the physics. In that case, use instead a **0-D PSR
(perfectly-stirred reactor) or plug-flow-reactor tabulation** (an ISAT-style
approach: pre-tabulate reactor solutions over your operating envelope,
interpolate at runtime) -- `cooling_pfr.py` in this kit already provides the
PFR building block for exactly this single-stream case; a premixed-flame
regime-1 zone would reuse that same idiom rather than `flamelet.py`.

Do not force the flamelet abstraction onto a premixed or single-stream
problem just because this kit ships it -- verify your actual regime-1
topology (two distinct feed streams that mix as they burn) before using
`Flamelet`.
