#!/usr/bin/env wolframscript

Remove["Global`*"]

(* Original equations *)
ueqn = \[Lambda][t]*uu[t]*(vv[t]-bb0)
veqn = vv[t]*(1-\[Kappa][t]*(uu[t]+vv[t])) - \[Lambda][t]*uu[t]*vv[t]

(* Lambda and kappa defintions *)
lt[t_] := (1/aa0)*(1-bb0*\[Kappa][t]) - \[Kappa][t]
l0 := (1/aa0)*(1-bb0*kk0) - kk0
\[Lambda][t_] := \[Alpha]*lt[t] + (1-\[Alpha])*l0
\[Kappa][t_] := kk0 + kk1*Cos[(1/n)*\[Omega]0*t]

(* Might use this later *)
(*\[Omega] = Sqrt[(b0/a0)*(b0*k0-1)^2+(1/4)*(3*b0*k0-4)] *)

(* rescaling away a0 *)
bb0 = aa0*b0
kk0 = (1/aa0)*k0
kk1 = (1/aa0)*k1
uu[t] = aa0*u[t]
vv[t] = aa0*v[t]
aa0 = 1

(* Simplify equations *)
simplifiedeqn = FullSimplify[{ueqn, veqn}]
Print["Simplified equations: "]
Print[simplifiedeqn[[1]]]
Print[simplifiedeqn[[2]]]

(* linearizing the equations *)
u[t] = 1 + \[Epsilon]*ul[t]
v[t] = b0 + \[Epsilon]*vl[t]
lineareqns = FullSimplify[Normal[Series[simplifiedeqn, {\[Epsilon], 0, 1}]] /. \[Epsilon] -> 1]
Print["Linearized equations: "]
Print[lineareqns[[1]]]
Print[lineareqns[[2]]]
