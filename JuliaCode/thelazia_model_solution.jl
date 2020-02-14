# ENV["PYTHON"] ="/home/sauldiazinfantevelasco/anaconda2/bin"
using Pkg
#
# Pkg.build("PyCall")
Pkg.add("DifferentialEquations")
Pkg.add("JSON")
Pkg.add("CSV")
Pkg.add("Revise")
Pkg.add("DataFrames")
Pkg.add("Gadfly")
Pkg.add("JLD")
Pkg.add("HDF5")
Pkg.add("PyPlot")
using DifferentialEquations
using JSON
using JLD, HDF5
#using PyCall
#pygui(true)
#using PyPlot
using Gadfly
using Revise
using DataFrames
#
function compute_r_zero(p)
    n_c = p[1]; lambda_f = p[2]; lambda_c = p[3];
    beta_c = p[4]; beta_c_tilde = p[5];
    beta_f = p[6]; beta_f_tilde = p[7];
    k_f = p[8]; k_c = p[9]; mu_f = p[10];
    mu_c = p[11]; n_f = lambda_f / mu_f
    #
    frac_1 = k_f / (mu_f + k_f)
    frac_2 = k_c / (mu_c + k_c)
    frac_3 = beta_c / mu_f
    frac_4 = (n_f / n_c) * (beta_f / mu_c)
    r_zero = (frac_1 * frac_2 * frac_3 * frac_4) ^ 0.25
    println("R0:\t", r_zero)
    return r_zero
end

function rhs(du, u, p, t)
    n_c = p[1]; lambda_f = p[2]; lambda_c = p[3];
    beta_c = p[4]; beta_c_tilde = p[4];
    beta_f = p[5];  beta_f_tilde = p[6];
    k_f = p[7]; k_c = p[8]; mu_f = p[9];
    mu_c = p[10]; rho = p[11]; theta = p[12];
    #
    w_f = 0.0; u_l = 0.0; v_h = 0.0;
    #
    s_f = u[1]; l_f = u[2]; i_f = u[3];
    s_c = u[4]; l_c = u[5];
    i_c_l = u[6]; i_c_h = u[7];
    t_c = u[8];

    new_s_f = (
                lambda_f
                    - beta_f / n_c * s_f * i_c_l
                    - beta_f_tilde / n_c * s_f * i_c_h
                    - mu_f * s_f
                    - w_f * s_f
            )
    new_l_f = (
                beta_f / n_c * s_f * i_c_l
                + beta_f_tilde / n_c * s_f * i_c_h
                - (k_f + mu_f) * l_f - w_f * l_f
            )
    new_i_f = k_f * l_f - mu_f * i_f - w_f * i_f
    #
    new_s_c = (lambda_c
                - beta_c / n_c * s_c * i_f
                - mu_c * s_c
                + rho * t_c
                + u_l * i_c_h
                )
    new_l_c = beta_c / n_c * s_c * i_f - (k_c + mu_c) * l_c
    new_i_c_l  = ( theta * k_c * l_c
                    - (mu_c + delta) * i_c_l
                    - beta_c_tilde / n_c * i_c_l * i_f
                    - u_l * i_c_h
                )
    new_i_c_h = ( (1.0 - theta ) * k_c * l_c
                    - mu_c * i_c_h
                    + beta_c_tilde / n_c * i_c_l * i_f
                    - v_h * i_c_h
    )
    new_t_c = - (rho + mu_c) * t_c + v_h * i_c_h
    #
    du[1] = new_s_f
    du[2] = new_l_f
    du[3] = new_i_f
    du[4] = new_s_c
    du[5] = new_l_c
    du[6] = new_i_c_l
    du[7] = new_i_c_h
    du[8] = new_t_c
end

function load_parameters(file_name)
    n_c = 1000
    lambda_f = 500
    beta_c = .001
    beta_c_tilde = .001
    beta_f = .001
    beta_f_tilde = .001
    k_f = 0.07142857142857142 # 1/14
    k_c = 0.02857142857142857
    mu_f = 0.022222222222222223
    mu_c = 0.000925925925925926
    lambda_c = n_c * mu_c
    rho = 0.025
    theta = 0.001
    parameterJSON = JSON.parse("""{
        "n_c": 1000,
        "lambda_f": 500,
        "lambda_c": 0.9259259259259259,
        "beta_c": 0.001,
        "beta_c_tilde": 0.001,
        "beta_f": 0.001,
        "beta_f_tilde": 0.001,
        "k_f": 0.07142857142857142,
        "k_c": 0.02857142857142857,
        "mu_f": 0.022222222222222223,
        "mu_c": 0.000925925925925926,
        "rho": 0.025,
        "theta": 0.7
        }""");

    parameterDict = Dict(#
        "n_c" => 1000.0, "lambda_f" => 500.0, "lambda_c" => 0.9259259259259259,
        "beta_c" => 0.001, "beta_c_tilde" => 0.001, "beta_f" => 0.001,
        "beta_f_tilde" => 0.001, "k_f" => 0.07142857142857142,
        "k_c" => 0.02857142857142857, "mu_f" => 0.022222222222222223,
        "mu_c" => 0.000925925925925926, "rho" => .025, "theta" => .001);
    string_data = JSON.json(parameterDict)
    #JSON.print(parameterJSON)
    f = open(file_name, "w") do j
        write(j, string_data)
    end
    # loading data
    global parameterDict = Dict()

    open(file_name, "r") do f
        global parameterDict
        dicttxt = readline(f)  # file information to string
        parameterDict = JSON.parse(dicttxt)
    end
    p = [   parameterDict["n_c"];
            parameterDict["lambda_f"];
            parameterDict["lambda_c"];
            parameterDict["beta_c"];
            parameterDict["beta_c_tilde"];
            parameterDict["beta_f"];
            parameterDict["beta_f_tilde"];
            parameterDict["k_f"];
            parameterDict["k_c"];
            parameterDict["mu_f"];
            parameterDict["mu_c"];
            parameterDict["rho"];
            parameterDict["theta"]
    ]
    return p;
end
path = string(pwd(), "/default_parameters.json")
p = load_parameters(path);
n_c = p[1]; lambda_f = p[2]; lambda_c = p[3];
beta_c = p[4]; beta_c_tilde = p[5];
beta_f = p[6]; beta_f_tilde = p[7];
k_f = p[8]; k_c = p[9]; mu_f = p[10]; mu_c = p[11];
rho = p[12]; theta = p[13];
u_zero = [ #
            lambda_f / mu_f - 100.0; 40.0; 60.0;
            n_c - 10.0; 2.0; 8.0; 0.0; 0.0
        ];
t_span = (0.0, 2000);
ode_problem = ODEProblem(rhs, u_zero, t_span, p)
sol = solve(ode_problem)
t = sol.t
solDataFrame=[sol.t, sol[1, :], sol[2, :], sol[3, :],
                sol[4, :], sol[5, :], sol[6, :], sol[7,:], sol[8,:]]

using CSV
solDataFrame = DataFrame(solDataFrame)
p1 = plot(solDataFrame, x=:x1, y=:x2,
    Guide.xlabel("time (days)"),
    Guide.ylabel("Suceptibles Flyes"))
p2 = plot(solDataFrame, x=:x1, y=:x3,
    Guide.xlabel("time (days)"),
    Guide.ylabel("Latent Flyes"))
p3 = plot(solDataFrame, x=:x1, y=:x4,
    Guide.xlabel("time (days)"),
    Guide.ylabel("Infected Flyes"))
#
p4 = plot(solDataFrame, x=:x1, y=:x5,
    Guide.xlabel("time (days)"),
    Guide.ylabel("Suceptibles Cows"))
p5 = plot(solDataFrame, x=:x1, y=:x6,
    Guide.xlabel("time (days)"),
    Guide.ylabel("Latent Cows"))
p6 = plot(solDataFrame, x=:x1, y=:x7,
    Guide.xlabel("time (days)"),
    Guide.ylabel("Low parasitic load Infected Cows"))
p7 = plot(solDataFrame, x=:x1, y=:x8,
        Guide.xlabel("time (days)"),
        Guide.ylabel("High parasitic load Infected Cows"))
p8 = plot(solDataFrame, x=:x1, y=:x8,
                Guide.xlabel("time (days)"),
                Guide.ylabel("Treated Cows"))
#
n_flyes = solDataFrame[:, 2 : 4]
n_cows = solDataFrame[:, 5 : 7]
#
time = solDataFrame[:, 1]
n_flyes = n_flyes.x2[:] + n_flyes.x3[:] + n_flyes.x4[:]
n_cows = n_cows.x5[:] + n_cows.x6[:] + n_cows.x7[:]
#
#
p7 = plot(x=time, y=n_flyes,
    Guide.xlabel("time (days)"),
    Guide.ylabel("Flyes CL"))
p8 = plot(x=time, y=n_cows,
    Guide.xlabel("time (days)"),
    Guide.ylabel("Cows CL"))
title(hstack(p1, p2, p3), "Flyes")
title(hstack(p4, p5, p6), "Cows")
plt0 = vstack(
        hstack(p1, p2, p3),
        #hstack(p4, p5, p6),
        hstack(p7)
        )
plt1 = vstack(
        # hstack(p1, p2, p3),
        hstack(p4, p5, p6, p7),
        hstack(p8)
        )
img0 = SVG("flyes_disease_dynamics.svg", 19cm, 11.74289cm)
img1 = SVG("cows_disease_dynamics.svg", 19cm, 11.74289cm)
draw(img0, plt0)
draw(img1, plt1)
