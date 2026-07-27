# EchTherm Shell-And-Tube Modeling Extraction

Generated from the PNG documentation pages and the added
`ravigururajan1996.pdf` paper in this folder on 2026-07-10.

This note is an engineering extraction, not a verbatim copy. It keeps the
correlations, variable definitions, validity ranges, and implementation-relevant
remarks needed to compare EchTherm behavior with the Combustor-HX model.

Source images:

- `image (17).png` to `image (8).png`: shell-side Bell-Delaware heat-transfer
  and pressure-drop documentation.
- `image (7).png` to `image.png`: internal corrugated/ribbed tube correlations.

Additional extracted source:

- `ravigururajan1996_extracted.md`: raw text extraction from
  `ravigururajan1996.pdf`, used to complete the Ravigururajan-Bergles friction
  and heat-transfer correlations that are partly cut off in `image (3).png`.

## 1. Shell-Side Heat Transfer

EchTherm uses a Bell-Delaware-style shell-side heat-transfer model. The real
shell-side coefficient is obtained from an ideal cross-flow coefficient corrected
by leakage, bypass, window-flow, and entrance/exit compartment factors:

```text
alpha_real = alpha_ideal * Jc * Jf * Jb * Js * Jr
```

where:

- `alpha_ideal`: ideal cross-flow heat-transfer coefficient across an equivalent
  tube bundle;
- `Jc`: correction for tubes in baffle windows;
- `Jf`: shell/baffle and tube/baffle leakage correction;
- `Jb`: bundle-shell bypass correction;
- `Js`: inlet/outlet baffle-spacing correction;
- `Jr`: adverse temperature-gradient correction, mostly important at low
  Reynolds number.

The ideal heat-transfer coefficient is expressed through the Colburn factor:

```text
alpha_ideal = j * Cp * Pr^(-2/3) * G * Ctherm
```

The documentation uses `G` as shell-side mass velocity:

```text
Re = G * Dt / mu = (m_dot / Sm) * Dt / mu
```

where:

- `Dt`: tube outside diameter;
- `Sm`: cross-flow area in the bundle;
- `mu`: dynamic viscosity;
- `m_dot`: shell-side mass flow.

The Colburn factor is:

```text
j = a1 * (1.33 / (Ltp / Dt))^a * Re^a2
a = a3 / (1 + 0.14 * Re^a4)
```

`Ltp` is the tube pitch. The coefficients depend on tube layout and Reynolds
number.

### 1.1 Ideal Heat-Transfer Coefficients

The extracted coefficient table is:

| Layout | a3 | a4 |
|---|---:|---:|
| Triangular 30 deg | 1.45 | 0.519 |
| Square 90 deg | 1.187 | 0.370 |
| Rotated square 45 deg | 1.93 | 0.5 |

| Reynolds range | Triangular 30 deg a1 | Triangular 30 deg a2 | Square 90 deg a1 | Square 90 deg a2 | Rotated 45 deg a1 | Rotated 45 deg a2 |
|---|---:|---:|---:|---:|---:|---:|
| Re < 10 | 1.4 | -0.667 | 0.97 | -0.667 | 1.55 | -0.667 |
| 10 <= Re < 100 | 1.36 | -0.657 | 0.9 | -0.631 | 0.498 | -0.656 |
| 100 <= Re < 1e3 | 0.593 | -0.477 | 0.408 | -0.46 | 0.73 | -0.5 |
| 1e3 <= Re < 1e4 | 0.321 | -0.388 | 0.107 | -0.266 | 0.37 | -0.396 |
| 1e4 <= Re < 1e5 | 0.321 | -0.388 | 0.37 | -0.395 | 0.37 | -0.396 |

### 1.2 Thermophysical Correction `Ctherm`

For liquids:

```text
Ctherm = (mu / mu_p)^0.14
```

where `mu_p` is viscosity evaluated at wall temperature. For cooling,
`Ctherm = 1`.

For gases during heating:

```text
Ctherm = ((T + 273) / (Tp + 273))^0.25
```

The documentation uses Celsius-form temperatures in this expression. In code,
the equivalent Kelvin expression should be used consistently:

```text
Ctherm = (T_bulk_K / T_wall_K)^0.25
```

### 1.3 Window Correction `Jc`

`Jc` accounts for the fact that flow in baffle windows is more parallel to the
tubes and has different velocity than the main cross-flow region.

```text
Jc = 0.55 + 0.72 * Fc
```

where `Fc` is the fraction of tubes in the cross-flow/baffled zone. The
documentation states that baffle cut/window opening should remain between about
15 percent and 45 percent of the shell diameter. It reports approximate `Jc`
behavior:

- no tubes in windows: about `Jc = 1.15`;
- small window: around `Jc = 1`;
- large window: down to about `Jc = 0.52`.

### 1.4 Leakage Correction `Jf`

`Jf` accounts for leakage through shell/baffle and tube/baffle clearances:

```text
Jf = 0.44 * (1 - Rsb) + [1 - 0.44 * (1 - Rsb)] * exp(-2.2 * Rlm)
```

with:

```text
Rsb = Ssb / (Ssb + Stb)
Rlm = (Ssb + Stb) / Sm
```

where:

- `Ssb`: shell-baffle leakage area;
- `Stb`: tube-baffle leakage area;
- `Sm`: cross-flow area in the bundle.

The documentation notes that shell-baffle leakage is thermally important because
it can bypass a large part of the heat-transfer process. Tube-baffle leakage is
usually less severe because heat transfer can still occur, but fouling/clogging
may matter. It says `Jf` should not be below about `0.6` and is classically
between `0.7` and `0.8`.

### 1.5 Bundle-Shell Bypass Correction `Jb`

`Jb` accounts for bundle-shell bypass flow:

```text
Jb = exp(-Cbht * Fsbp * (1 - (2 * Rss)^(1/3)))
```

and:

```text
Jb = 1 for Rss > 0.5
```

with:

```text
Fsbp = Sb / Sm
Rss = Nss / Ntcc
```

where:

- `Sb`: bypass area;
- `Sm`: cross-flow area in the bundle;
- `Nss`: number of sealing-strip pairs;
- `Ntcc`: number of tube rows in the cross-flow zone.

Thermal factor `Cbht`:

```text
Cbht = 1.35  for laminar flow
Cbht = 1.25  for turbulent flow, Re >= 100
```

The documentation remarks that typical `Jb` values depend strongly on exchanger
head type and diametral bundle-shell clearance. Minimum expected `Jb` values are
around `0.7`; lower values suggest adding sealing strips to force more flow
through the bundle.

### 1.6 Entrance/Exit Compartment Correction `Js`

For equal inlet and outlet compartment lengths:

```text
Js = ((Nb - 1) + 2 * Lad^(1 - Nad)) / ((Nb - 1) + 2 * Lad)
```

with:

```text
Lad = Lio / Lbc
```

where:

- `Nb`: number of baffles;
- `Lio`: inlet/outlet baffle compartment length;
- `Lbc`: central baffle spacing;
- `Nad`: exponent or slope parameter from the correction curve.

The documentation notes that inlet and outlet compartments are generally longer
than central compartments, commonly by less than a factor of about 2.

The supplied screenshots do not explicitly state the numerical value of `Nad`.
The maintained Combustor-HX Bell-Delaware implementation uses the standard
turbulent heat-transfer exponent:

```text
Nad = 0.6
```

For unequal inlet and outlet spacings this is implemented in the more general
form:

```text
Js = [(Nb - 1) + (Lsi/Lbc)^(1 - Nad) + (Lso/Lbc)^(1 - Nad)]
     / [(Nb - 1) + (Lsi/Lbc) + (Lso/Lbc)]
```

For equal inlet/outlet spacings, `Lsi = Lso = Lio`, this reduces to the EchTherm
form shown above.

### 1.7 Laminar Adverse-Temperature-Gradient Correction `Jr`

The image set shows `Jr` in the global Bell-Delaware heat-transfer product but
does not include the page defining it. Public search did not locate a freely
available primary equation page. The candidate equation below is the standard
Bell-Delaware/Serth form already used in the local Combustor-HX implementation;
it should be treated as the comparison target unless the missing EchTherm page
is later recovered.

```text
Jr = 1                                                for Re_s >= 100

Jr_20 = (10 / Nc)^0.18
Jr = Jr_20                                           for Re_s <= 20

Jr = Jr_20 + (1 - Jr_20) * (Re_s - 20) / 80          for 20 < Re_s < 100
```

where:

- `Re_s`: shell-side Reynolds number based on tube outside diameter and
  shell-side maximum cross-flow mass velocity;
- `Nc`: total number of effective tube rows crossed by the shell-side stream.

Combustor-HX currently evaluates:

```text
Nc = (N_tcc + N_tcw) * (Nb + 1)
```

This is conservative in the sense that `Jr` is exactly `1` for the normal
turbulent shell-side regime and only matters for low-Re shell-side operation.

## 2. Shell-Side Pressure Drop

EchTherm decomposes shell-side pressure drop into:

```text
DeltaP_total = DeltaP_c + DeltaP_e + DeltaP_w
```

where:

- `DeltaP_c`: cross-flow pressure drop in central baffle compartments;
- `DeltaP_e`: entrance/exit compartment pressure drop;
- `DeltaP_w`: window pressure drop.

The documentation computes all three from two base hydraulic configurations:

- ideal cross-flow pressure drop, `DeltaP_bi`;
- window pressure drop, `DeltaP_wi`.

Then it applies the same correction philosophy as Bell-Delaware:

```text
DeltaP_c = DeltaP_bi * (Nb - 1) * Rb * Rl
DeltaP_e = DeltaP_bi * (1 + Ntcw / Ntcc) * Rb * Rs
DeltaP_w = DeltaP_wi * Nb * Rl
```

where:

- `Rb`: bundle-shell bypass correction for pressure drop;
- `Rl`: leakage correction for pressure drop;
- `Rs`: inlet/outlet spacing correction for pressure drop;
- `Ntcw`: number of tube rows in a window;
- `Ntcc`: total number of tube rows in the cross-flow zone.

### 2.1 Ideal Cross-Flow Pressure Drop

```text
DeltaP_bi = 2 * f * Ntcc * G^2 / (rho * Ctherm)
```

The friction factor is:

```text
f = b1 * (1.33 / (Ltp / Dt))^b * Re^b2
b = b3 / (1 + 0.14 * Re^b4)
```

`b1` and `b2` depend on Reynolds number and layout. `b3` and `b4` depend only on
layout.

### 2.2 Ideal Pressure-Drop Coefficients

The extracted coefficient table is:

| Layout | b3 | b4 |
|---|---:|---:|
| Triangular 30 deg | 7.0 | 0.5 |
| Square 90 deg | 6.3 | 0.378 |
| Rotated square 45 deg | 6.59 | 0.52 |

| Reynolds range | Triangular 30 deg b1 | Triangular 30 deg b2 | Square 90 deg b1 | Square 90 deg b2 | Rotated 45 deg b1 | Rotated 45 deg b2 |
|---|---:|---:|---:|---:|---:|---:|
| Re < 10 | 48 | -1 | 35 | -1 | 32 | -1 |
| 10 <= Re < 100 | 45.1 | -0.973 | 32.1 | -0.963 | 26.2 | -0.913 |
| 100 <= Re < 1e3 | 4.57 | -0.476 | 6.09 | -0.602 | 3.5 | -0.476 |
| 1e3 <= Re < 1e4 | 0.486 | -0.152 | 0.0815 | 0.022 | 0.333 | -0.136 |
| 1e4 <= Re < 1e5 | 0.372 | -0.123 | 0.391 | -0.148 | 0.303 | -0.126 |

### 2.3 Window Pressure Drop

For laminar flow, `Re < 100`:

```text
DeltaP_wi = 26 * (Gw * mu / rho) * (Ntcw / (Ltp - Dt) + Lbc / Dw^2) + Gw^2 / rho
```

For turbulent flow, `Re > 100`:

```text
DeltaP_wi = (2 + 0.6 * Ntcw) * Gw^2 / (2 * rho)
```

where:

- `Ntcw`: number of tube rows in a baffle window considered as cross-flow rows;
- `Ltp`: tube pitch;
- `Dt`: tube diameter;
- `Lbc`: central baffle spacing;
- `Dw`: hydraulic diameter of the window;
- `Gw`: window mass velocity.

The documentation defines `Gw` as total mass flow divided by the geometric mean
of the cross-flow and window flow areas:

```text
Gw = m_dot / sqrt(Sm * Sw)
```

where `Sw` is the window flow area.

### 2.4 Leakage Pressure-Drop Correction `Rl`

```text
Rl = exp[-1.33 * (1 + Rsb) * Rlm^P]
P = -0.15 * (1 + Rsb) + 0.8
```

with:

```text
Rsb = Ssb / (Ssb + Stb)
Rlm = (Ssb + Stb) / Sm
```

### 2.5 Bundle-Shell Bypass Pressure-Drop Correction `Rb`

```text
Rb = exp[-Cbp * Fsbp * (1 - (2 * Rss)^(1/3))]
Rb = 1 for Rss > 0.5
```

with:

```text
Fsbp = Sb / Sm
Rss = Nss / Ntcc
```

`Cbp` depends on flow regime:

```text
Cbp = 4.5  for laminar flow
Cbp = 3.7  for turbulent flow, Re > 100
```

The documentation states that this factor also accounts for different inlet and
outlet compartment lengths relative to central compartments.

### 2.6 Entrance/Exit Pressure-Drop Correction `Rs`

For equal inlet/outlet compartment lengths:

```text
Rs = 2 * Lad^(n - 2)
Lad = Lio / Lbc
```

where `n` is the slope of the friction curve:

```text
n = 1.0  for laminar flow
n = 0.2  for turbulent flow
```

## 3. Internal Corrugated / Ribbed Tube Correlations

The EchTherm documentation also contains internal tube-side correlations for
single-phase flow in corrugated or ribbed tubes. These are relevant to the
EchTherm geometry-page option "Grooved Tube".

The documentation distinguishes:

- corrugated tubes: corrugation depth comparable to, or only moderately smaller
  than, tube diameter;
- ribbed/rainured tubes: corrugation depth much smaller than tube diameter, with
  many corrugations.

The geometric parameters are:

- `D`: tube diameter measured at the bottom of corrugation;
- `theta`: helix/corrugation angle relative to the tube axis, between 0 deg and
  80 deg;
- `e`: corrugation/rib height, typically 0.2 mm to 2 mm depending on internal
  tube diameter;
- `P`: axial pitch between two successive ribs/corrugations;
- `N`: number of helical starts, from 1 to 10 for corrugated tubes and possibly
  above 10 for ribbed tubes.

The documentation states:

- Nakayama is mainly for circular/sinusoidal corrugations;
- Ravigururajan and Bergles covers more general corrugation types;
- Webb is dedicated to ribbed tubes.

### 3.1 Nakayama Method

Nakayama introduces two functions, `R` and `G`, to evaluate roughness effects on
near-wall heat/momentum transfer.

For heat transfer:

```text
G = (lambda / (8 * St) - 1) / sqrt(lambda / 8) + R
```

For friction:

```text
R = 1 / sqrt(lambda / 8) + 2.5 * ln(2e / D) + 3.75
```

Stanton number:

```text
St = h / (rho * U * Cp)
```

Pressure drop:

```text
DeltaP = lambda * (L / D) * rho * U^2 / 2
```

The documentation classifies flow near corrugations by angle:

| Corrugation angle | Regime description |
|---|---|
| theta > 60 deg | transverse flow near the wall |
| 45 deg < theta < 60 deg | mixed transition between transverse and swirling/turbulent |
| theta < 45 deg | swirling/turbulent flow, mostly parallel to corrugations |

It defines a roughness Reynolds number:

```text
e_plus = (e / D) * Re * sqrt(lambda / 8)
```

For transverse flow, `theta > 60 deg`:

```text
R = 4.5 + 5.63e-4 * (P/e)^2.59 * ln(e_plus)
G = 4.75 * e_plus^0.28 * Pr^0.57
```

Validity notes:

```text
11 < P/e < 14  for R
11 < P/e < 57  for G
```

For swirling/turbulent flow, `theta < 45 deg`:

```text
R = 5.02 * e_plus^0.15 * (theta / 45)^(-0.16) * (P * sin(theta) / e)^0.1
G = 4.75 * e_plus^0.28 * Pr^0.57
```

Validity:

```text
11 < P/e < 57
```

For transition flow, `45 deg < theta < 60 deg`:

```text
R = 5.14 * e_plus^0.12 * (theta / 45)^(-0.8) * (P * sin(theta) / e)^0.1
```

The corresponding `G` expression in the screenshot is:

```text
G = 4.9 * e_plus^0.37 * Pr^0.57
```

with validity:

```text
15 < P/e < 19  for R
14 < P/e < 18  for G
15 < e_plus < 200
```

The documentation describes an iterative algorithm:

1. Define `P/e`, `theta`, and `D`.
2. Compute Reynolds and Prandtl numbers.
3. Initialize `lambda_i = 0.02`.
4. Compute `e_plus` from `lambda_i`.
5. Compute `R` for the selected corrugation-angle regime.
6. Compute a new `lambda` from the friction relation and iterate until
   convergence.
7. Compute `G`.
8. Compute Stanton number and finally the convective coefficient `h`.

### 3.2 Ravigururajan And Bergles

This correlation is valid for water and air in the following ranges:

```text
5000 <= Re <= 250000
0.66 <= Pr <= 37.6
0.01 <= e/D <= 0.2
0.1 <= P/D <= 7
0.3 <= alpha/90 <= 1
```

The enhanced-tube Nusselt number is given relative to smooth-tube Nusselt:

```text
Nu_a / Nu_s =
{1 + [2.64 * Re^0.036 * (e/D)^0.212 * (P/D)^(-0.21)
      * (alpha/90)^0.29 * Pr^0.024]^7}^{1/7}
```

where:

- `Nu_a`: corrugated/enhanced tube Nusselt number;
- `Nu_s`: smooth tube Nusselt number.

The smooth-tube baseline `Nu_s` is given by the Petukhov correlation:

```text
Nu_s = ((f_s / 2) * Re * Pr)
       / (1 + 12.7 * sqrt(f_s / 2) * (Pr^(2/3) - 1))
```

Here `f_s` is the smooth-tube Fanning friction coefficient. This is important
because the Ravigururajan-Bergles Nusselt relation is an enhancement ratio, not
a standalone absolute Nusselt model.

The smooth-tube Fanning friction factor is the Filonenko expression:

```text
f_s = [1.58 * ln(Re) - 3.28]^(-2)
```

The full final enhanced-tube friction correlation is visible in the
Ravigururajan-Bergles paper and mostly visible in the EchTherm screenshot. The
paper labels it Eq. (4); the EchTherm page labels the same form Eq. (14):

```text
f_a / f_s =
{1 + [
    29.1 * Re^(0.67 - 0.06 * P/D - 0.49 * alpha/90)
    * (e/D)^(1.37 - 0.157 * P/D)
    * (P/D)^(-1.66e-6 * Re - 0.33 * alpha/90)
    * (alpha/90)^(4.59 + 4.11e-6 * Re - 0.15 * P/D)
    * (1 + 2.94 * sin(beta) / n)
  ]^(15/16)}^(16/15)
```

where:

- `f_a`: augmented/enhanced-tube Fanning friction factor;
- `f_s`: smooth-tube Fanning friction factor;
- `alpha`: helix angle of the rib/corrugation;
- `beta`: contact angle or profile angle used in the shape function;
- `n`: number of sharp corners facing the flow; the paper states `n = 2` for
  triangular and rectangular ribs and `n -> infinity` for smoother profiles.

The paper first defines the empirical shape-function idea as:

```text
Sh = (1 + const / n) * sin(beta)
```

then states that the fitted constant is `2.94`. The final printed paper
equation and the EchTherm screenshot use the compact factor:

```text
1 + 2.94 * sin(beta) / n
```

Use that compact factor for implementation. The earlier ambiguity in this
extraction came from reading the prose definition and the printed final equation
side by side; the final equation is the implementation-relevant one.

The documentation notes that, unlike Nakayama, this method directly computes an
enhanced-tube Fanning friction factor and pressure drop as:

```text
DeltaP = 2 * f_a * (L / D) * rho * U^2
```

The paper reports that the friction correlation predicts 96 percent of the
database within +/-50 percent and 77 percent within +/-20 percent; the
heat-transfer correlation predicts 99 percent within +/-50 percent and
69 percent within +/-20 percent. It also states that the heat-transfer
correlation intentionally has no shape function because profile shape did not
show a marked heat-transfer influence in the database.

### 3.3 Webb Correlation

The Webb correlation is for ribbed tubes and is valid primarily in:

```text
20000 <= Re <= 80000
5.08 <= Pr <= 6.29  for water
25 <= N <= 45
25 deg <= alpha <= 45 deg
D = 15.54 mm
0.33 mm <= e <= 0.55 mm
```

The Colburn factor is:

```text
j = 0.00933 * Re^(-0.181) * N^0.285 * (e/D)^0.323 * theta^0.505
```

The friction factor is:

```text
f = 0.108 * Re^(-0.283) * N^0.221 * (e/D)^0.785 * theta^0.78
```

The documentation then relates pressure drop to Fanning friction:

```text
lambda = 4f = 2 * DeltaP * D / (rho * L * U^2)
```

Additional validity notes:

```text
25 deg < theta < 45 deg
18 < N < 45
0.021 < e/D < 0.0356
5.08 < Pr < 6.29
20000 < Re < 80000 for the Colburn j factor
12000 < Re < 80000 for the friction factor
```

## 4. References Extracted From The Screenshots

The screenshots provide nine references in total: four for the Bell-Delaware
shell-side model and five for the corrugated/ribbed tube models.

### 4.1 Shell-Side / Bell-Delaware References

1. Bell K. J., 1963, *Final report of the cooperative research program on shell
   and tube heat exchangers*, Bulletin 5, University of Delaware Engineering
   Experiment Station, New York.

   DOI/source status: no DOI found in Crossref; no OpenLibrary catalog record
   found from title/author queries. Current traceable source is the EchTherm
   screenshot reference itself.

2. Bell K. J., 1981, Delaware method for shellside design, in Kakac S.,
   Bergles A. E., and Mayinger F. eds., *Heat Exchangers: Thermal Hydraulic
   Fundamentals and Design*, pp. 581-618.

   DOI/source status: no standalone chapter DOI found in Crossref. Book source
   record: OpenLibrary `Heat Exchangers`, work
   https://openlibrary.org/works/OL8636692W .

3. Bell K. J., 1988, Delaware method for shellside design, in Shah R. K.,
   Subbarao E. C., and Mashelkar R. A. eds., *Heat Transfer Equipment Design*,
   pp. 145-166.

   DOI/source status: no standalone chapter DOI found in Crossref. Book source
   record: OpenLibrary `Heat transfer equipment design`, work
   https://openlibrary.org/works/OL19493326W .

4. Taborek J., 1969, shell and tube heat exchanger design - extension of the
   method to other shell, baffle, and tube bundle geometries, *Heat Exchanger
   Design Handbook*, vol. 3, section 3.3.11.

   DOI/source: https://doi.org/10.1615/hedhme.a.000257 .

### 4.2 Corrugated / Ribbed Tube References

5. Nakayama W., 1983, Spiral Ribbing to enhance single phase heat transfer inside
   tube, *ASME/JSME Conference*.

   DOI/source status: no DOI found in Crossref for the cited conference item.
   Current traceable source is the EchTherm screenshot reference itself.

6. Li H. M., 1982, Investigation on tube side flow visualization, *7th
   International Heat Transfer Conference*.

   DOI/source: https://doi.org/10.1615/IHTC7.1500 .

7. Withers J. G., 1980, Tube side heat transfer and pressure drop for tubes having
   internal spiral rib with turbulent/transitional flow of single phase fluids,
   *Heat Transfer Engineering*, vol. 2, no. 1-2.

   DOI/source: Crossref returns two same-title `Heat Transfer Engineering`
   records, likely the two-part article in vol. 2, no. 1-2:
   https://doi.org/10.1080/01457638008962755 and
   https://doi.org/10.1080/01457638008962750 .

8. Ravigururajan T. S. and Bergles A. E., 1996, Development and Verification of
   General Correlations for Pressure Drop and Heat Transfer in Single-Phase
   Turbulent Flow in Enhanced Tubes, *Experimental Thermal and Fluid Science*,
   13:55-70, New York.

   DOI/source: https://doi.org/10.1016/0894-1777(96)00014-3 .

9. Webb R. L., Narayanamurthy R., and Thors P., 2000, Heat transfer and friction
   characteristics of internal helical-rib roughness, *ASME Journal of Heat
   Transfer*, 122(1):134-142.

   DOI/source: https://doi.org/10.1115/1.521444 .

## 5. Immediate Comparison Targets For Combustor-HX

This section is intentionally limited to comparison targets, not conclusions.

For the shell-side model, compare whether Combustor-HX currently has:

- the same ideal cross-flow Colburn `j` form and layout/Re tables;
- the same pressure-drop `f` form and layout/Re tables;
- explicit `Jc`, `Jf`, `Jb`, `Js`, and `Jr` corrections;
- explicit `Rl`, `Rb`, and `Rs` pressure-drop corrections;
- window pressure drop with separate laminar/turbulent formulas;
- thermophysical correction `Ctherm` for gases and liquids.

For the grooved tube model, compare whether Combustor-HX currently has:

- Nakayama angle-regime logic and iterative friction/Stanton solution;
- Ravigururajan-Bergles Nusselt enhancement and friction enhancement;
- Webb ribbed-tube `j` and `f` options;
- validity warnings when the EchTherm geometry lies outside each correlation
  range.

## 6. Screenshot Coverage Audit

The currently available screenshots were rechecked after adding the
Ravigururajan-Bergles paper.

Covered explicitly from the screenshots:

- shell-side global heat-transfer product
  `alpha_real = alpha_ideal * Jc * Jf * Jb * Js * Jr`;
- ideal cross-flow heat transfer, Colburn factor, layout/Re coefficient tables,
  and `Ctherm`;
- `Jc`, `Jf`, `Jb`, and `Js`;
- shell-side pressure-drop decomposition into central cross-flow,
  entrance/exit, and window components;
- ideal cross-flow pressure drop and pressure-drop coefficient tables;
- laminar and turbulent window pressure-drop formulas;
- hydraulic `Rl`, `Rb`, and `Rs`;
- Nakayama method, including `R`, `G`, `St`, `DeltaP`, `e_plus`, the three
  angle regimes, and the iterative solution workflow;
- Ravigururajan-Bergles validity range, Nusselt enhancement, Petukhov baseline,
  Filonenko friction, and final friction enhancement;
- Webb Colburn/friction correlations and validity ranges;
- all reference lists visible in the supplied screenshots.

Recovered from added sources or standard implementation context:

- Ravigururajan-Bergles final friction equation is completed from
  `ravigururajan1996.pdf`.
- The implementation-relevant Ravigururajan-Bergles shape factor is resolved as
  `1 + 2.94 * sin(beta) / n`.
- Heat-transfer `Js` exponent `Nad = 0.6` and the unequal inlet/outlet spacing
  generalization are added from the maintained Combustor-HX Bell-Delaware
  implementation and standard Bell-Delaware convention.
- `Jr` low-Re correction is added from the maintained Combustor-HX
  Bell-Delaware implementation and standard Bell-Delaware/Serth convention.

Still not fully available in the screenshots:

- the detailed expression for `Jr` is not present in the image set. The global
  product includes `Jr`, but the folder does not contain the page that defines
  its formula or validity handling. The formula now included in Section 1.7 is
  therefore a standard-method completion, not a screenshot extraction.
