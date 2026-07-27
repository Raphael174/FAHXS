""" 
@ author : Raphael Aubry

Spray and evaporation characteristics of ethanol and gasoline direct
injection in non-evaporating, transition and flash-boiling conditions

https://sci-hub.se/10.1016/j.enconman.2015.10.081
"""

import matplotlib.pyplot as plt
import numpy as np 

#%% data from paper on gasoline vapour pressures vs temperature

{
"x": [275.7891,280.3662,284.6275,287.7841,292.3611,297.4116,301.3573,305.4609,310.3535,314.7727,318.7184,322.9798,326.9255,331.3447,335.1326,338.447,342.5505,347.4432,352.3359,356.9129,361.6477,366.2247,370.0126,373.4848,377.4306,382.6389,387.5316,391.7929,395.7386,399.3687],
"y": [9.5949,11.7271,14.9254,15.9915,20.2559,24.5203,29.8507,34.1151,40.5117,45.8422,51.1727,59.7015,66.0981,75.693,87.42,94.8827,105.5437,119.403,136.4606,153.5181,172.7079,191.8977,213.2196,229.2111,249.467,280.3838,312.3667,343.2836,368.8699,398.7207]
}



def fit_exponential_from_data(x, y):
    """
    Fit y ≈ a * exp(b * (x - x0))  using least squares on ln(y).
    Returns params dict and an evaluator function f(xq).

    Notes:
      - y must be > 0 for all points (log is taken).
      - x0 centers x to improve numerical conditioning.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size != y.size or x.size < 2:
        raise ValueError("x and y must be same-length arrays (len >= 2).")
    if np.any(y <= 0.0):
        raise ValueError("All y values must be positive for exponential fit.")

    x0 = float(x.mean())                     # center x to reduce correlation
    X = x - x0
    ln_y = np.log(y)

    # Linear least squares: ln(y) = ln(a) + b * (x - x0)
    # Using polyfit for readability; equivalent to lstsq on [1, X].
    b, ln_a = np.polyfit(X, ln_y, 1)
    a = float(np.exp(ln_a))

    # Build an evaluator
    def f(xq):
        xx = np.asarray(xq, dtype=float)
        yhat = a * np.exp(b * (xx - x0))
        return yhat if np.ndim(xq) else float(yhat)

    # Simple quality metric (R^2) for sanity check
    y_fit = f(x)
    ss_res = float(np.sum((y - y_fit)**2))
    ss_tot = float(np.sum((y - y.mean())**2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    params = {"a": a, "b": float(b), "x0": x0, "r2": r2}
    return params, f

x = [275.7891,280.3662,284.6275,287.7841,292.3611,297.4116,301.3573,305.4609,310.3535,314.7727,318.7184,322.9798,326.9255,331.3447,335.1326,338.447,342.5505,347.4432,352.3359,356.9129,361.6477,366.2247,370.0126,373.4848,377.4306,382.6389,387.5316,391.7929,395.7386,399.3687]
y = [9.5949,11.7271,14.9254,15.9915,20.2559,24.5203,29.8507,34.1151,40.5117,45.8422,51.1727,59.7015,66.0981,75.693,87.42,94.8827,105.5437,119.403,136.4606,153.5181,172.7079,191.8977,213.2196,229.2111,249.467,280.3838,312.3667,343.2836,368.8699,398.7207]

params, y_exp = fit_exponential_from_data(x, y)
print("Fit params:", params)
print("y at x=360:", y_exp(360.0))   # scalar
xs = np.linspace(min(x), max(x), 5)
print("batch:", y_exp(xs))           # vector


x_ = np.linspace(275, 400, 300, endpoint=True)
y_=y_exp(x_)

plt.scatter(x, y, label="original data")
plt.plot(x_, y_, label="curve fit")
plt.legend()
plt.show()