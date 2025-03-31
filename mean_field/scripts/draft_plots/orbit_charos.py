# %%
import numpy as np
from matplotlib import pyplot as plt
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
orbit02 = np.loadtxt('../../output/time_series/alpha_0.2_orbit.dat')
orbit0 = np.loadtxt('../../output/time_series/alpha_0_orbit.dat')

a02 = orbit02[40000:, 1]
a0 = orbit0[40000:, 1]
b02 = orbit02[40000:, 2]
b0 = orbit0[40000:, 2]
t02 = orbit02[40000:, 0]
t0 = orbit0[40000:, 0]

us = 1
k0 = 0.25
k1 = 0.24
n = 0.7
alpha02 = 0.2
alpha0 = 0
afp02 = afp(us, k0, k1, n, alpha02, t02)
afp0 = afp(us, k0, k1, n, alpha0, t0)
# %%
plt.rcParams.update({'font.size': 20})  # You can change 14 to any size you prefer

plt.figure()
plt.plot(a02, b02, c='black', label='Orbit')
plt.plot(afp02, np.ones_like(afp02), c = 'red', label='Fixed point')
plt.xlabel('a')
plt.ylabel('b', rotation=0, labelpad=20)
plt.legend()
plt.savefig('../../plots/draft_chaos_bifurcation/alpha_0.2_orbit.pdf', dpi=100, bbox_inches='tight')
plt.show()

# %%
plt.figure()
plt.plot(a0, b0, c='black', label='Orbit')
plt.plot(afp0, np.ones_like(afp0), c = 'red', label='Fixed point')
plt.xlabel('a')
plt.ylabel('b', rotation=0, labelpad=20)
plt.legend()
plt.savefig('../../plots/draft_chaos_bifurcation/alpha_0_orbit.pdf', dpi=100, bbox_inches='tight')
plt.show()
