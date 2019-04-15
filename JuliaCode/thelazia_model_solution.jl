ENV["PYTHON"] ="/home/sauldiazinfantevelasco/anaconda2/bin"
using Pkg
#Pkg.build("PyCall")
Pkg.add("DifferentialEquations")
Pkg.add("JSON")
Pkg.add("PyCall")
# Pkg.add("Gadfly")
# Pkg.add("JLD")
# Pkg.add("HDF5")
Pkg.add("PyPlot")
using DifferentialEquations
using JSON
using JLD, HDF5
pygui(true)
using PyPlot


function rhs(du, u, p, t)
    n_c = p[1]; lambda_f = p[2]; lambda_c = p[3];
    beta_c = p[4]; beta_f = p[5]; k_f = p[6];
    k_c = p[7]; mu_f = p[8]; mu_c = p[9];
    #
    s_f = u[1]; l_f = u[2]; i_f = u[3];
    s_c = u[4]; l_c = u[5]; i_c = u[6];
    new_s_f = lambda_f - beta_f / n_c * s_f * i_c - mu_f * s_f
    new_l_f = beta_f / n_c * s_f * i_c - (k_f + mu_f) * i_f
    new_i_f = k_f * i_f - mu_f * i_f
    #
    new_s_c = lambda_f - beta_c / n_c * s_c * i_f - mu_c * s_c
    new_l_c = beta_c / n_c * s_c * i_f - (k_c + mu_c) * i_c
    new_i_c = k_c * i_c - mu_c * i_c

    du[1] = new_s_f
    du[2] = new_l_f
    du[3] = new_i_f
    du[4] = new_s_c
    du[5] = new_l_c
    du[6] = new_i_c
end

function load_parameters(file_name)
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

    parameterDict = Dict(#
        "n_c" => 1000.0, "lambda_f" => 500.0, "lambda_c" => 200.0,
        "beta_c" => 0.001, "beta_f" => 0.001, "k_f" => 0.07142857142857142,
        "k_c" => 0.02857142857142857, "mu_f" => 0.022222222222222223,
        "mu_c" => 0.000925925925925926);
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
            parameterDict["beta_f"];
            parameterDict["k_f"];
            parameterDict["k_c"];
            parameterDict["mu_f"];
            parameterDict["mu_c"]
    ]
    return p;
end
path = string(pwd(),"/default_parameters.json")
p = load_parameters(path);
n_c = p[1]; lambda_f = p[2]; lambda_c = p[3];
beta_c = p[4]; beta_f = p[5]; k_f = p[6];
k_c = p[7]; mu_f = p[8]; mu_c = p[9];
u_zero = [n_c - 10.0; 2.0; 8.0; lambda_f / mu_f - 100.0; 40.0; 60.0];
t_span = (0.0, 100);
ode_problem = ODEProblem(rhs, u_zero, t_span, p)
sol = solve(ode_problem)
plot(sol.t, sol[1,:])
