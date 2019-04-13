
using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("JSON")
Pkg.add("JLD")
using DifferentialEquations
using JSON
using JLD


function rhs(du, u, p, t)
    s_f = u[1]
    l_f = u[2]
    i_f = u[3]
    s_c = u[4]
    l_c = u[5]
    i_c = u[6]
    new_s_f = lambda_f - beta_f / n_c * s_f * i_c - mu_f * s_f
    new_l_f = beta_f / n_c * s_f * i_c - (k_f + mu_f) * i_f
    new_i_f = k_f * i_f - mu_f * i_f
    #
    new_s_c = lambda_f - beta_c / n_c * s_c * i_f - mu_c * s_c
    new_l_c = beta_c / n_c * s_c * i_f - (k_c + mu_c) * i_c
    new_i_c = k_c * i_c - mu_c * i_c

    du = [s_f; l_f; i_f; s_c; l_c; i_c]
end

n_c = 1000
lambda_f = 500
lambda_c = 200
beta_c = .001
beta_f = .001
k_f = 0.07142857142857142 # 1/14
k_c = 0.02857142857142857
mu_f = 0.022222222222222223
mu_c = 0.000925925925925926

parameterJSON = JSON.parse("""{
    "n_c": 1000,
    "lambda_f": 500,
    "lambda_c": 200,
    "beta_c": 0.001,
    "beta_f": 0.001,
    "k_f": 0.07142857142857142,
    "k_c": 0.02857142857142857,
    "mu_f": 0.022222222222222223,
    "mu_c": 0.000925925925925926
    }""");
JSON.print(stdout, parameterJSON)
open("default_parameters.json", "w") do j
    write(j, parameterJSON)
end
