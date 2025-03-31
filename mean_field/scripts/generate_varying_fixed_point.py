# %%
import numpy as np
import matplotlib.pyplot as plt
# %%
# Functions
def o0(us, k0):
    a1 = 4*us
    a1_9 = 1-k0
    a2 = 4*a1_9**2-4*us*k0*a1_9-us*k0**2
    o0 = np.sqrt(a2/a1)
    return o0

def kappa(us, k0, k1, o, t):
    return k0 + k1*np.cos(o*t)

def lamba(us, k0, k1, n, alpha, t):
    lamba0 = (1-(1+us)*k0)/us
    o = o0(us, k0)/n
    lamba1 = (1-(1+us)*kappa(us, k0, k1, o, t))/us
    return (1-alpha)*lamba0 + alpha*lamba1

def afp(us, k0, k1, n, alpha, t):
    o = o0(us, k0)/n
    return (1-kappa(us, k0, k1, o, t))/(lamba(us, k0, k1, n, alpha, t)+kappa(us, k0, k1, o, t))
# %%
# Parameters
us = 1
k0 = 0.25
k1 = 0.24
n = 0.7
# %%
plt.figure()
t = np.linspace(0, 50, 1000)
alpha = 0.2
afp_values = afp(us, k0, k1, n, alpha, t)
plt.plot(afp_values, np.ones_like(afp_values), label='alpha = 0.2')
plt.show()
