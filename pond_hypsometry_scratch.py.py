# Scratchwork for sea level pond parameterization

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.optimize import root, curve_fit

rhosi = 940.0
rhofresh = 1000.0
rhow = 1026.0

a_p_star = 0.25
hi = np.linspace(0.1, 5, 100)
hi = np.array([0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0])

# Calculations for linear hypsometry
r = hi*(rhow - rhosi)/(rhofresh*a_p_star**2 - 2*rhow*a_p_star + rhow)
hp = r * a_p_star
hsp = hi - r + 2*hp
hocn = (hi*rhosi + a_p_star*hp*rhofresh)/rhow

# Calculations for logistic hypsometry
def e(a, L, k, a_0, y):
    """Logistic hypsometric function for elevation above ice base
    
    Parameters
    ----------
    a : float
        Fractional area of category
    L : float
        range of elevations (max - min)
    k : float
        steepness of hypsometric curve
    a_0 : float
        a value of midpoint of logistic curve
    y : float
        vertical offset of elevation curve
    
    Returns
    -------
    float
        value of hypsometric curve at a
    """

    return L/(1 + np.exp(-k*(a - a_0))) + y

def E(a, L, k, a_0, y):
    """Integral of hypsometric function
    
    Parameters
    ----------
    a : float
        Fractional area of category
    L : float
        range of elevations (max - min)
    k : float
        steepness of hypsometric curve
    a_0 : float
        a value of midpoint of logistic curve
    y : float
        vertical offset of elevation curve
    
    Returns
    -------
    float
        value of integral of hypsometric curve at a
    """

    return L*np.log(np.exp(-k*(a - a_0)) + 1)/k + (L + y)*a - L*np.log(np.exp(k*a_0) + 1)/k

def hocn(hi, Vp, rhosi, rhop, rhow):
    """Compute ocean height above mean base of the ice assuming rigid isostacy
    
    Parameters
    ----------
    hi : float
        Mean ice thickness (m)
    Vp : float
        Pond volume per unit area of category (m)
    rhosi : float
        Sea ice density (kg/m3)
    rhop : float
        Pond water density (kg/m3)
    rhow : float
        Sea water density (kg/m3)
    
    Returns
    -------
    float
        Ocean height above ice base (m)
    """

    return (hi*rhosi + Vp*rhop)/rhow

def Vp_a(a, L, k, a_0, y):
    """Compute pond volume given a pond area and logistic parameters
    
    Parameters
    ----------
    a : float
        pond fraction per unit category area (m)
    L : float
        range of elevations (max - min)
    k : float
        steepness of hypsometric curve
    a_0 : float
        a value of midpoint of logistic curve
    y : float
        vertical offset of elevation curve

    Returns
    -------
    float
        pond volume per unit category area (m)
    """

    return a*e(a, L, k, a_0, y) - E(a, L, k, a_0, y)

def dVpda(a, L, k, a_0, y):
    """Compute derivative of pond volume with respect to a
    
    Parameters
    ----------
    a : float
        pond fraction per unit category area (m)
    L : float
        range of elevations (max - min)
    k : float
        steepness of hypsometric curve
    a_0 : float
        a value of midpoint of logistic curve
    y : float
        vertical offset of elevation curve

    Returns
    -------
    float
        derivative of pond volume per unit category area w.r.t. a
    """

    return a*k*(e(a,L,k,a_0,y) - y)*(1 - (e(a,L,k,a_0,y) - y)/L)

def ap(Vp, L, k, a_0, y, hi, ap_0):
    """Compute pond area given a pond volume and logistic parameters
    
    Parameters
    ----------
    Vp : float
        Pond volume per unit category area (m)
    L : float
        range of elevations (max - min)
    k : float
        steepness of hypsometric curve
    a_0 : float
        a value of midpoint of logistic curve
    y : float
        vertical offset of elevation curve
    hi : float
        mean category ice thickness (m)
    ap_0 : float
        initial guess for pond fraction
    
    Returns
    -------
    float
        pond area fraction of category
    """

    tol = 0.001

    # Check if pond fills 100% of area
    if (Vp >= Vp_a(1.0, L, k, a_0, y)):
        return 1.0
    else: # Newton's method
        ap = ap_0
        V = Vp_a(ap, L, k, a_0, y)
        while (np.abs(V-Vp) > tol):
            ap = ap - (Vp_a(ap, L, k, a_0, y) - Vp)/dVpda(a, L, k, a_0, y)
    return ap

# Given target sea level pond area fraction and depth ap' hp', k, and hi
# find the parameters a_0, L, and y that satisfy the following constraints
# E(1) - hi = 0 # mean elevation equals ice thickness
# ap'*e(ap') - E(ap') - ap'*hp' = 0 ! target pond dimensions are consistent with hypsometric curve
# e(ap') - hocn = 0 ! pond surfaces at sea level
# Implement these constraints in a function that we can use with fsolve
def obj_fun(vec, *constants):
    """Objective function
    
    Parameters
    ----------
    vec : array (3,)
        Array of parameters to optimize [L, a_0, y]
    constants : tuple
        tuple of constant values to unpack

    Returns
    -------
    array (3,)
        Right hand side values of the three constraints    
    """

    # Unpack constants
    ap_sl, hp_sl, k, hi, rhosi, rhop, rhow = constants

    out = np.zeros(3)

    # mean elevation constraint
    out[0] = E(1, vec[0], k, vec[1], vec[2]) - hi
    # target pond area/depth matching hypsometry constraint
    out[1] = (ap_sl*e(ap_sl, vec[0], k, vec[1], vec[2]) 
              - E(ap_sl, vec[0], k, vec[1], vec[2])
              - ap_sl * hp_sl)
    # pond surface at sea level constraint
    out[2] = (e(ap_sl, vec[0], k, vec[1], vec[2]) 
              - hocn(hi, ap_sl*hp_sl, rhosi, rhop, rhow))
    return out

# Test out function
# targets
ap_sl = 0.28
hp_sl = 0.14
# params
k = 15
hi = 10 #2 #0.2 #0.5 #1.0#2.0
rhosi = 940.0
rhofresh = 1000.0
rhop = rhofresh
rhow = 1026.0

args = (ap_sl, hp_sl, k, hi, rhosi, rhop, rhow)

# initial guess for [L, a_0, y]
L_g = 0.1* hi #0.02# 1.0#0.125 #0.25#0.5
# ap_sl*y + (1 - ap_sl)*(y+L_g) = hi
# y (ap_sl + 1 - ap_sl) + L_g - L_g * ap_sl = hi
# CHECK !y = hi - L_g*(1 - ap_sl)
guess = np.array([L_g, ap_sl, hi - L_g*(1 - ap_sl)])

sol = root(obj_fun, guess, args=args)
# For this to produce good results we, need k > 5ish
a = np.linspace(0, 1, 50)
e_guess = e(a, guess[0], k, guess[1], guess[2])
e_res = e(a, sol.x[0], k, sol.x[1], sol.x[2])
hocn_res = hocn(hi, ap_sl*hp_sl, rhosi, rhop, rhow)
plt.plot(a, e_res)
plt.plot(a, e_guess, '--')
plt.hlines(y=hocn_res, xmin=0, xmax=1, linestyles=':')
plt.show()

# Compute coefficients for a range of hi
num = 50
hi_min = 0.5
hi_max = 30
his = np.geomspace(hi_min, hi_max, num=num)
Ls = np.zeros(num)
a_0s = np.zeros(num)
ys = np.zeros(num)
for i in np.arange(num):
    L_g = 0.1 * his[i]
    guess = np.array([L_g, ap_sl, his[i] - L_g*(1 - ap_sl)])
    args = (ap_sl, hp_sl, k, his[i], rhosi, rhop, rhow)
    sol = root(obj_fun, guess, args=args)
    Ls[i] = sol.x[0]
    a_0s[i] = sol.x[1]
    ys[i] = sol.x[2]

f, axs = plt.subplots(1, 3, figsize=(15, 5))
axs[0].plot(his, Ls, '.')
axs[1].plot(his, a_0s, '.')
axs[2].plot(his, ys, '.')

# Plot results
a = np.linspace(0, 1)
f, axs = plt.subplots(1, 2, figsize= (12, 7))
for i in np.arange(num):
    e_res = e(a, Ls[i], k, a_0s[i], ys[i])
    axs[0].plot(a, e_res, c=cm.viridis((his[i]-hi_min)/(hi_max - hi_min)))
    hocn_res = hocn(his[i], ap_sl*hp_sl, rhosi, rhop, rhow)
    axs[1].plot(a, e_res-hocn_res, c=cm.viridis((his[i]-hi_min)/(hi_max - hi_min)))

axs[1].hlines(0, xmin=0, xmax=1, colors='k', linestyles="--")
norm = Normalize(vmin=hi_min, vmax=hi_max)

# Play around with curves to fit these functions
# Fit Ls, polynomial
def f_L(hi, *params):
    hi_0, c02, c0, c1, c2 = params 
    return c02*(hi - hi_0)**(1/2) + c0 + c1*(hi - hi_0) + c2*(hi - hi_0)**(2)
hi_0 = 0.49
c02 = 0.5
c0 = 0.32
c1 = 1/6
c2 = 0.01
p0 = (hi_0, c02, c0, c1, c2)
opt_L, cov_L = curve_fit(f_L, his, Ls, p0=p0)
plt.plot(his, Ls, '+')
#plt.plot(his, f_L(his, *p0))
plt.plot(his, f_L(his, *opt_L))
print(np.sqrt(((f_L(his, *opt_L)-Ls)**2).sum()))
# Fit a_0, polynomial
def f_a_0(hi, *params):
    hi_0, c0, c1, c2, c3 = params
    return c0 + c1*(hi - hi_0) + c2*(hi - hi_0)**(1/2) + c3*(hi - hi_0)**(1/3)
c1 = .002
hi_0 = 0.49
c2 = .05
c3 = 0.1
c0 = 0.2
p0 = (hi_0, c0, c1, c2, c3)
opt_a_0, cov_a_0 = curve_fit(f_a_0, his, a_0s, p0=p0)
plt.plot(his, a_0s, '+')
#plt.plot(his, f_a_0(his, *p0))
plt.plot(his, f_a_0(his, *opt_a_0))
print(np.sqrt(((f_a_0(his, *opt_a_0)-a_0s)**2).sum()))
# Fit y, polynomial
def f_y(hi, *params):
    hi_0, c03, c02, c0, c1 = params 
    return c03*(hi - hi_0)**(1/3) + c02*(hi - hi_0)**(1/2) + c0 + c1*(hi - hi_0)
hi_0 = 0.49
c03 = 1
c02 = 1
c0 = 0.5
c1 = 1
p0 = (hi_0, c03, c02, c0, c1)
opt_y, cov_y = curve_fit(f_y, his, ys, p0=p0)
plt.plot(his, ys, '+')
plt.plot(his, f_y(his, *opt_y))
print(np.sqrt(((f_y(his, *opt_y)-ys)**2).sum()))
