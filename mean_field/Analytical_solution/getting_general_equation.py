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
# %%
# Find the eigenvectors and the p_matrix
zero_order_eqs = [eqs[0].subs({a: 0, k1: 0}), eqs[1].subs({a: 0, k1: 0})]
jac = sp.Matrix([[zero_order_eqs[0].rhs.coeff(u), zero_order_eqs[0].rhs.coeff(v)], [zero_order_eqs[1].rhs.coeff(u), zero_order_eqs[1].rhs.coeff(v)]])
x = sp.Matrix([u, v])
zero_order_eq = sp.Eq(x.diff(t), jac*x)
eigen_sys = jac.eigenvects()
eigen_value = eigen_sys[0][0]
# p_matrix diagonalizes the system
p_matrix = sp.Matrix.hstack(eigen_sys[0][2][0], eigen_sys[1][2][0])
p_matrix = p_matrix.applyfunc(sp.simplify)
p_matrix_inv = p_matrix.inv()
p_matrix_inv = p_matrix_inv.applyfunc(sp.simplify)
# j_diag is the daionalized_jacobian
j_diag = p_matrix_inv*jac*p_matrix
j_diag = j_diag.applyfunc(sp.simplify)
# %%
# Seperate dampening and frequency for zero order and write everything in terms of zeta and omega_2
damp = eigen_value.args[1]
freq_sq = -eigen_value.args[0]**2
z, o02 = sp.symbols('zeta omega_2')
temp_eq_1 = sp.Eq(damp, z)
temp_eq_2 = sp.Eq(freq_sq, o02)
sol = sp.solve([temp_eq_1, temp_eq_2], (vs, k0))
vs_in_terms_of_z_o02 = sol[0][0]
k0_in_terms_of_z_o02 = sol[0][1]
eigen_value_in_terms_of_z_o02 = eigen_value.subs({vs: vs_in_terms_of_z_o02, k0: k0_in_terms_of_z_o02}).simplify()
p_matrix_in_terms_of_z_o02 = p_matrix.subs({vs: vs_in_terms_of_z_o02, k0: k0_in_terms_of_z_o02})
p_matrix_in_terms_of_z_o02 = p_matrix_in_terms_of_z_o02.applyfunc(sp.simplify)
p_matrix_inv_in_terms_of_z_o02 = p_matrix_inv.subs({vs: vs_in_terms_of_z_o02, k0: k0_in_terms_of_z_o02})
p_matrix_inv_in_terms_of_z_o02 = p_matrix_inv_in_terms_of_z_o02.applyfunc(sp.simplify)
zero_order_eq_in_terms_of_z_o02 = zero_order_eq.subs({vs: vs_in_terms_of_z_o02, k0: k0_in_terms_of_z_o02})
zero_order_eq_in_terms_of_z_o02 = zero_order_eq_in_terms_of_z_o02.simplify()
jac_in_terms_of_z_o02 = jac.subs({vs: vs_in_terms_of_z_o02, k0: k0_in_terms_of_z_o02})
jac_in_terms_of_z_o02 = jac_in_terms_of_z_o02.applyfunc(sp.simplify)
j_diag_in_terms_of_z_o02 = j_diag.subs({vs: vs_in_terms_of_z_o02, k0: k0_in_terms_of_z_o02})
j_diag_in_terms_of_z_o02 = j_diag_in_terms_of_z_o02.applyfunc(sp.simplify)
# %%
# Transform coordinates to eigenvectors
x_trans = p_matrix_inv_in_terms_of_z_o02*x
x_trans_diff = x_trans.diff(t).subs({u.diff(t): zero_order_eq_in_terms_of_z_o02.rhs[0], v.diff(t): zero_order_eq_in_terms_of_z_o02.rhs[1]})
x_trans_diff = x_trans_diff.applyfunc(sp.simplify)
# NOTE: All of the above is correct and checked
x_in_terms_of_x_trans = p_matrix_in_terms_of_z_o02*x
# %%
x_trans_diff_incorrect = x_trans_diff.subs({u: x_in_terms_of_x_trans[0], v: x_in_terms_of_x_trans[1]})
x_trans_diff_incorrect = x_trans_diff_incorrect.applyfunc(sp.simplify)
