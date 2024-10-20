# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---
# %%
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import sympy as sp
from scipy.signal import find_peaks
# %%
# Numerical solution
def LV_zero(t, y, m, l):
    a, b = y
    dadt = -m*a + l*a*b
    dbdt = b*(1-b/k) - l*a*b
    return [dadt, dbdt]
# %%
m = 0.01
l = 0.1
k = 1
# fixed point is (1/l, m/l)
y0_deviation = [1,1]
y0 = [(1/l)*(1-(m/(l*k)))+y0_deviation[0], (m/l)+y0_deviation[1]]
t_eval = np.linspace(0, 500, 5000)
# %%
solution = solve_ivp(LV_zero, [0, 500], y0, args=(m, l), t_eval=t_eval)
t = solution.t 
a = solution.y[0]
b = solution.y[1]
# %%
plt.figure(figsize=(10, 6))
plt.plot(t, a, label='a')
plt.plot(t, b, label='b')
plt.legend()
plt.show()
# %%
plt.plot(a, b, label='phase space plot')
plt.xlabel('a')
plt.ylabel('b')
plt.legend()
plt.show()
# %%
peaks, _ = find_peaks(a)
periods = np.diff(t[peaks])
average_period = np.mean(periods)
frequency = 1 / average_period
print(frequency*2*np.pi, average_period)
# %%
if len(peaks) > 1:
    midpoints = []
    avg_a_periods = []
    avg_b_periods = []
    
    for i in range(1, len(peaks)):
        start_idx = peaks[i - 1]
        end_idx = peaks[i]
        midpoints.append((t[start_idx] + t[end_idx]) / 2)
        avg_a_periods.append(np.mean(a[start_idx:end_idx]))
        avg_b_periods.append(np.mean(b[start_idx:end_idx]))

    # Plot the results
    plt.figure(figsize=(12, 8))
    plt.plot(midpoints, avg_a_periods, label='Average a over each period', marker='o')
    plt.plot(midpoints, avg_b_periods, label='Average b over each period', marker='x')
    plt.xlabel('Time (midpoint of periods)')
    plt.ylabel('Average values')
    plt.title('Average of a and b over each period as a function of time')
    plt.legend()
    plt.show()
else:
    print("Not enough peaks to compute averages over periods.")
# %%
fft_a = np.fft.fft(a)
fft_b = np.fft.fft(b)
fft_a[0] = 0
fft_b[0] = 0
sample_spacing = t[1] - t[0]
frequencies = np.fft.fftfreq(a.size, d=sample_spacing)
amplitude_a = np.abs(fft_a) / len(a)
amplitude_b = np.abs(fft_b) / len(b)
# %%
plt.figure(figsize=(12, 6))
plt.plot(frequencies[:len(frequencies)//2], amplitude_a[:len(frequencies)//2], label='Amplitude of a')
plt.plot(frequencies[:len(frequencies)//2], amplitude_b[:len(frequencies)//2], label='Amplitude of b')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('Fourier Transform of a and b')
plt.legend()
plt.grid(True)
plt.show()



# %%
# Analytical solution
# TODO: perturbation expansion no longer works, might have to reparametrize the equations in terms of the fixed point
t_s, l_s, m_s, k_s = sp.symbols('t lambda mu k', positive=True, real=True)
a_s, b_s = sp.symbols('a b', cls=sp.Function)
a_s = a_s(t_s)
b_s = b_s(t_s)
a_eq = sp.Eq(a_s.diff(t_s), -m_s*a_s + l_s*a_s*b_s)
b_eq = sp.Eq(b_s.diff(t_s), b_s * (1-(b_s/k_s)) - l_s*a_s*b_s)
# %%
# Equation for fluctuations around fixed point
fixed_point = sp.solve([a_eq.rhs, b_eq.rhs], (a_s, b_s))
fixed_point = fixed_point[2]
subs = {a_s: a_s + fixed_point[0], b_s: b_s + fixed_point[1]}
a_eq = a_eq.subs(subs).simplify()
b_eq = b_eq.subs(subs).simplify()
# %%
# Expand the equation to n-th order in lambda
n = 1
a_terms = [sp.symbols(f'a{i}', cls=sp.Function)(t_s) for i in range(n + 1)]
b_terms = [sp.symbols(f'b{i}', cls=sp.Function)(t_s) for i in range(n + 1)]
a_expanded = sum(l_s**i * a_terms[i] for i in range(n + 1))
b_expanded = sum(l_s**i * b_terms[i] for i in range(n + 1))
subs = {a_s: a_expanded, b_s: b_expanded}
a_rhs_expanded = a_eq.rhs.subs(subs).series(l_s, n=n+1).removeO()
b_rhs_expanded = b_eq.rhs.subs(subs).series(l_s, n=n+1).removeO()
a_rhs_coeffs = [a_rhs_expanded.coeff(l_s, i) for i in range(n + 1)]
b_rhs_coeffs = [b_rhs_expanded.coeff(l_s, i) for i in range(n + 1)]
# %%
a0, b0 = sp.symbols('a_i b_i', real=True, positive=True)
solution_s = []
for i in range(n+1):
    # Solve the equation to n-th order in lambda
    a_eq_order_n = sp.Eq(a_terms[i].diff(t_s), a_rhs_coeffs[i])
    b_eq_order_n = sp.Eq(b_terms[i].diff(t_s), b_rhs_coeffs[i])
    ics = {a_terms[i].subs(t_s, 0): a0 if i == 0 else 0, b_terms[i].subs(t_s, 0): b0 if i == 0 else 0}
    solution_s.append(sp.dsolve([a_eq_order_n, b_eq_order_n], ics=ics))
    solution_s[i][0] = solution_s[i][0].simplify()
    solution_s[i][1] = solution_s[i][1].simplify()
    if i < n:
        a_rhs_coeffs[i+1] = a_rhs_coeffs[i+1].subs({a_terms[i]: solution_s[i][0].rhs, b_terms[i]: solution_s[i][1].rhs})
        b_rhs_coeffs[i+1] = b_rhs_coeffs[i+1].subs({b_terms[i]: solution_s[i][1].rhs, a_terms[i]: solution_s[i][0].rhs})
        a_rhs_coeffs[i+1] = a_rhs_coeffs[i+1].simplify()
        b_rhs_coeffs[i+1] = b_rhs_coeffs[i+1].simplify()
# %%
a_final = sp.Integer(0)
b_final = sp.Integer(0)
for i in range(n+1):
    a_final += l**i * solution_s[i][0].rhs
    b_final += l**i * solution_s[i][1].rhs
a_solution = sp.lambdify(t_s, a_final.subs({l_s: l, m_s: m, a0: y0_deviation[0], b0: y0_deviation[1]}), 'numpy')
b_solution = sp.lambdify(t_s, b_final.subs({l_s: l, m_s: m, a0: y0_deviation[0], b0: y0_deviation[1]}), 'numpy')
a_analytical = (1/l) + a_solution(t)
b_analytical = (m/l) + b_solution(t)
plt.figure(figsize=(10, 6))
plt.plot(t, a_analytical, label='a analytical')
plt.plot(t, a, label='a numerical')
plt.legend()
plt.show()

