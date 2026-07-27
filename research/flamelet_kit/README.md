# flamelet_kit

A standalone, portable extraction of the RIF (Representative Interactive
Flamelet) finite-rate-chemistry methodology, retargeted from a rocket-engine
ignition/extinction problem to a **steady or steadily-changing** two-regime
heat exchanger: a non-premixed mixing/flame zone feeding a downstream
already-mixed, cooling hot-gas channel (e.g. diesel/O2 injector -> heat
exchanger channel).

No ignition/extinction machinery, no NN surrogates, no telemetry. Depends on
only `numpy`, `scipy`, and `cantera`. Zero dependency on the source
`hybrid_rocket` package this was extracted from.

## The two regimes

```
  Regime 1: MIXING / FLAME ZONE          Regime 2: COOLING CHANNEL
  ---------------------------            --------------------------
  Fuel stream  \                          One already-mixed hot-gas
                >-- non-premixed flame --> stream flows down a channel,
  Oxidizer stream /   (flamelet.py)        cooling via wall heat loss
                                           (cooling_pfr.py)

  Two feed streams, a conserved           No second stream, no mixture
  scalar Z is meaningful ->               fraction -> plug-flow reactor
  flamelet-in-mixture-fraction-space      (0-D reactor marched along x)
```

`Flamelet` (regime 1) solves the flame; `CoolingPFR` (regime 2) takes the
flame's burnt-gas state as its inlet and marches finite-rate chemistry with
wall heat loss down the channel. See `METHODOLOGY.md` for the physics and
`ADAPTATION_GUIDE.md` for how each regime couples to the other and how to
retarget both to your own mechanism/fuel/geometry.

## Install

```
pip install numpy scipy cantera
```

## 30-second quickstart

```python
import numpy as np, cantera as ct
from flamelet_kit import Flamelet, SteadyCache, CoolingPFR

# --- Regime 1: the flame -------------------------------------------------
gas = ct.Solution("gri30.yaml")
gas.TPX = 300., 101325., "CH4:1";            Y_fuel = np.array(gas.Y)
gas.TPX = 300., 101325., "O2:1, N2:3.76";    Y_ox   = np.array(gas.Y)

fl = Flamelet("gri30.yaml", n_z=65)
fl.init_mixing(T_ox=300., Y_ox=Y_ox, T_fuel=300., Y_fuel=Y_fuel, p=101325.)
for _ in range(200):
    fl.step(dt=2e-5, p=101325., T_ox=300., Y_ox=Y_ox,
            T_fuel=300., Y_fuel=Y_fuel, chi_st=3.0)
print(fl.T_max, fl.Z_st)

# --- Cache: skip re-advancing when conditions haven't moved --------------
cache = SteadyCache()
if cache.is_hit(p_key=101325., chi_key=3.0, T_fuel=300.):
    pass  # reuse fl as-is
else:
    fl.step(...)  # re-advance, then cache.record_advance(...)

# --- Regime 2: the cooling channel, fed by regime 1's burnt state -------
T_burnt = fl.T_at_Z(fl.Z_st)
Y_burnt = np.array([fl.Y_at_Z(fl.Z_st, sp) for sp in fl.species_names])
pfr = CoolingPFR("gri30.yaml")
result = pfr.march(T_in=T_burnt, Y_in=Y_burnt, p=101325., mdot=0.02,
                    length=1.0, n_steps=40, diameter=0.02,
                    h_conv=200., T_wall=500.)
print(result["T"][-1], result["Q_wall_total_W"])
```

Run the full runnable demos:

```
python flamelet_kit/example_run.py       # regime 1: flamelet
python flamelet_kit/example_cooling.py   # regime 2: cooling PFR
python -m pytest flamelet_kit/tests/ -q  # test suite
```

## File map

| File | What |
|---|---|
| `flamelet.py` | `Flamelet` class (regime 1: non-premixed flame in mixture-fraction space) + grid/Bilger-Z/chi-profile/CN-diffusion helpers |
| `steady_cache.py` | `SteadyCache`: tolerance-gated "is the stored solution still valid" reuse decision -- the main cost lever, for either regime |
| `cooling_pfr.py` | `CoolingPFR` class (regime 2: single-stream plug-flow reactor with wall heat loss) |
| `flamelet_bank.py` | Optional multi-condition manager pattern; **most users need only ONE `Flamelet`, skip this** (see its module docstring) |
| `example_run.py` | Runnable regime-1 demo (GRI-30 CH4/air, bundled with Cantera) |
| `example_cooling.py` | Runnable regime-2 demo, feeding an adiabatic-flame stand-in state into `CoolingPFR` |
| `tests/test_flamelet_kit.py` | pytest suite for both regimes |
| `METHODOLOGY.md` | The physics: mixture-fraction space, Peters flamelet equations, chi(Z), Strang split, PFR + wall heat loss, steady-cache rationale |
| `ADAPTATION_GUIDE.md` | How to retarget: mechanism/fuel swap (incl. diesel surrogate note), chi_st estimation routes, cache tolerance tuning, what was removed and why, applicability caveat |
| `REPRODUCE_SPEC.md` | From-scratch algorithmic spec (no code) for both modules, precise enough to re-implement without reading the Python |

## Applicability (read before use)

`Flamelet` (regime 1) is for **non-premixed** (two-feed-stream,
diffusion-controlled) reacting flow where a conserved mixture fraction Z is
physically meaningful -- e.g. a diesel/O2 injector where fuel and oxidizer
arrive separately and mix as they burn. If your regime-1 zone is instead
**premixed** or single-stream, the flamelet-in-Z model is the wrong reduced
model; use a 0-D PSR / plug-flow-reactor tabulation (the same idiom
`cooling_pfr.py` already provides) instead. See `ADAPTATION_GUIDE.md` for
the full discussion.
