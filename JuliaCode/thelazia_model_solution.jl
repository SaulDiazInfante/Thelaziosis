# ENV["PYTHON"] ="/home/sauldiazinfantevelasco/anaconda2/bin"
using Pkg
#
# Pkg.build("PyCall")
Pkg.add("Cairo")
Pkg.add("DifferentialEquations")
Pkg.add("JSON")
Pkg.add("CSV")
Pkg.add("Revise")
Pkg.add("DataFrames")
Pkg.add("Gadfly")
Pkg.add("JLD")
Pkg.add("HDF5")
Pkg.add("PyPlot")
Pkg.add("Fontconfig")
using DifferentialEquations
using JSON
using JLD, HDF5
#using PyCall
#pygui(true)
#using PyPlot
using Cairo
using Fontconfig
using Gadfly
using Revise
using DataFrames
include("thelazia_model.jl")
#

path = string(pwd(), "/default_parameters.json")
p = load_parameters(path);
n_c_inf = p[1];
lambda_f = p[2];
lambda_c = p[3];
beta_c = p[4];
beta_c_tilde = p[5];
beta_f = p[6];
beta_f_tilde = p[7];
k_f = p[8];
k_c = p[9];
mu_f = p[10];
mu_c = p[11];
rho = p[12];
theta = p[13];
u_zero = [ #
            lambda_f / mu_f - 100.0; 40.0; 60.0;
            n_c_inf - 50.0; 2.0; 8.0; 20.0; 0.0
        ];
t_span = (0.0, 10000);
ode_problem = ODEProblem(uncontrollted_rhs, u_zero, t_span, p)
sol = solve(ode_problem)
t = sol.t
solDataFrame=[sol.t, sol[1, :], sol[2, :], sol[3, :],
                sol[4, :], sol[5, :], sol[6, :],
                 sol[7, :], sol[8, :]]
using CSV
solDataFrame = DataFrame(solDataFrame)
p1 = plot(solDataFrame, x=:x1, y=:x2,
    Guide.xlabel("time (days)"),
    Guide.ylabel("S_fs"), Geom.path)
p2 = plot(solDataFrame, x=:x1, y=:x3,
    Guide.xlabel("time (days)"),
    Guide.ylabel("L_f"), Geom.path)
p3 = plot(solDataFrame, x=:x1, y=:x4,
    Guide.xlabel("time (days)"),
    Guide.ylabel("I_f"), Geom.path)
#
p4 = plot(solDataFrame, x=:x1, y=:x5,
    Guide.xlabel("time (days)"),
    Guide.ylabel("S_c"), Geom.path)
p5 = plot(solDataFrame, x=:x1, y=:x6,
    Guide.xlabel("time (days)"),
    Guide.ylabel("L_c"), Geom.path)
p6 = plot(solDataFrame, x=:x1, y=:x7,
    Guide.xlabel("time (days)"),
    Guide.ylabel("I_cl"), Geom.path)
p7 = plot(solDataFrame, x=:x1, y=:x8,
        Guide.xlabel("time (days)"),
        Guide.ylabel("I_ch"), Geom.path)
p8 = plot(solDataFrame, x=:x1, y=:x9,
            Guide.xlabel("time (days)"),
            Guide.ylabel("T_c"), Geom.path)
#
n_flyes = solDataFrame[:, 2 : 4]
n_cows = solDataFrame[:, 5 : 9]
#
time = solDataFrame[:, 1]
n_flyes = n_flyes.x2[:] + n_flyes.x3[:] + n_flyes.x4[:]
n_cows = n_cows.x5[:] + n_cows.x6[:] + n_cows.x7[:] +
        + n_cows.x8[:] + n_cows.x9[:]
#
#
p9 = plot(x=time, y=n_flyes,
    Guide.xlabel("time (days)"),
    Guide.ylabel("Flyes CL"), Geom.path)
p10 = plot(x=time, y=n_cows,
    Guide.xlabel("time (days)"),
    Guide.ylabel("Cows CL"), Geom.path)
title(hstack(p1, p2, p3), "Flyes")
title(hstack(p4, p5, p6, p7, p8), "Cows")
plt0 = vstack(
        hstack(p1, p2, p3),
        #hstack(p4, p5, p6),
        hstack(p9)
        )
plt1 = vstack(
        # hstack(p1, p2, p3),
        hstack(p4, p5),
        hstack(p6, p7),
        hstack(p8, p10)
        )
img0 = PDF("flyes_disease_dynamics.pdf", 19cm, 11.74289cm)
img1 = PDF("cows_disease_dynamics.pdf", 19cm, 11.74289cm)
draw(img0, plt0)
draw(img1, plt1)
