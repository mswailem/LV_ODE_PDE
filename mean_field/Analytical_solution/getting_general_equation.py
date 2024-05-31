# %%
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
# %%
# Define the variables
k0, k1, vs, a, t, l0, m0, o0, n = sp.symbols('k_0 k_1 v_s alpha t lambda_0 mu_0 omega_0 n')
u = sp.Function('u')(t)
v = sp.Function('v')(t)
k = sp.Function('k')(t)
l = sp.Function('lambda')(t)
m = sp.Function('mu')(t)
la = sp.Function('lambda')(t)
ma = sp.Function('mu')(t)
l0 = 1-vs * k0 - k0
m0 = vs * l0
# o0 = sp.sqrt(vs*(vs*k0-1)**2 + (1/4)*vs*k0*(3*vs*k0-4))
k = k0 + k1 * sp.cos((o0/n)*t)
la = 1-vs*k-k
ma = vs*la
l = a*la + (1-a)*l0
m = a*ma + (1-a)*m0
# %%
# Define the general equation
eqs = [sp.Eq(u.diff(t), l*u*v-m*u), sp.Eq(v.diff(t), v*(1-k*(u+v)) - l*u*v)]
eqs[0] = eqs[0].simplify()
eqs[1] = eqs[1].simplify()
# %%
# Linearizing
eqs[0] = eqs[0].subs({u:1+u, v:vs+v}).simplify()
eqs[1] = eqs[1].subs({u:1+u, v:vs+v}).simplify()
eqs[0] = eqs[0].expand().subs({u**2: 0, u*v:0, v**2:0}).simplify()
eqs[1] = eqs[1].expand().subs({u**2: 0, u*v:0, v**2:0}).simplify()
