using Pkg
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
using CSV
using Gadfly
using DataFrames
# using Fontconfig
include("thelazia_model.jl")
# Simulation parameters
n_max = 10000; t_f = 2000.0
t_span = range(0.0, t_f, length=n_max)
x_dim = 8
u_dim = 3
u_control = zeros((n_max, u_dim))
x_path = zeros((n_max, x_dim))
psi_adjoint = zeros((n_max, x_dim))
h = t_span[2]

path = string(pwd(), "/default_parameters.json")
p = load_parameters(path);

# Forward sweep
function runge_kutta_forward(u)
    s_f_zero = p[13];
    l_f_zero = p[14];
    i_f_zero = p[15];
    s_c_zero = p[16];
    l_c_zero = p[17];
    i_c_l_zero = p[18];
    i_c_h_zero = p[19];
    t_c_zero = p[20];
    x_0 = [
                s_f_zero; l_f_zero; i_f_zero;
                s_c_zero; l_c_zero;
                i_c_l_zero; i_c_h_zero;
                t_c_zero
        ];
    x = zeros(n_max, x_dim)
    x[1, :] = x_0
    #
    for i = 1:(n_max - 1)
        x_i = x[i, :]
        u_next = u[i + 1, :]
        u_i = u[i, :]
        u_mean = 0.5 * (u_i + u_next)
        k_1 = rhs_f(x_i, u_mean)
        k_2 = rhs_f(x_i + 0.5 * k_1, u_mean)
        k_3 = rhs_f(x_i + 0.5 * k_2, u_mean)
        k_4 = rhs_f(x_i + 0.5 * k_3, u_mean)

        x_next = h / 6.0 * (k_1 + 2 * k_2 + 2 * k_3 + k_4)
        x[i + 1, :] = x_next
    end
    return x
end

function runge_kutta_backward(x)
    psi_s_f_zero = 0.0;
    psi_l_f_zero = 0.0;
    psi_i_f_zero = 0.0;
    psi_s_c_zero = 0.0;
    psi_l_c_zero = 0.0;
    psi_i_c_l_zero = 0.0;
    psi_i_c_h_zero = 0.0;
    psi_t_c_zero = 0.0;
    psi__0 = [
        psi_s_f_zero;
        psi_l_f_zero;
        psi_i_f_zero;
        psi_s_c_zero;
        psi_l_c_zero;
        psi_i_c_l_zero;
        psi_i_c_h_zero;
        psi_t_c_zero;
        ];
    psi = zeros(n_max, x_dim)
    psi[n_max, :] = psi_0
    #
    for i in n_max: -1: 2
        psi_i = psi[i, :]
        u_i = u_control[i, :]
        x_i = x[i, :]
        x_previous = x[i-1, :]
        x_mean = 0.5 * (x_i + x_previous)
        k_1 = rhs_adjoints(x_mean, u_i, psi_i)
        k_2 = rhs_f(x_i + 0.5 * k_1, u_mean)
        k_3 = rhs_f(x_i + 0.5 * k_2, u_mean)
        k_4 = rhs_f(x_i + 0.5 * k_3, u_mean)

        x_next = h / 6.0 * (k_1 + 2 * k_2 + 2 * k_3 + k_4)
        x[i + 1, :] = x_next
    end
    return x
end

function forward_plot()
    x_path = runge_kutta_forward(u_control)
    xDataFrame=[t_span, x_path[:, 1], x_path[:, 2], x_path[:, 3],
                    x_path[:, 4], x_path[:, 5], x_path[:, 6],
                    x_path[:, 7], x_path[:, 8]]
    solDataFrame = DataFrame(xDataFrame)
    p1 = plot(solDataFrame, x=:x1, y=:x2,
        Guide.xlabel("time (days)"),
        Guide.ylabel("S_f"), Geom.line)
    p2 = plot(solDataFrame, x=:x1, y=:x3,
        Guide.xlabel("time (days)"),
        Guide.ylabel("L_f"), Geom.line)
    p3 = plot(solDataFrame, x=:x1, y=:x4,
        Guide.xlabel("time (days)"),
        Guide.ylabel("I_f"), Geom.line)
    #
    p4 = plot(solDataFrame, x=:x1, y=:x5,
        Guide.xlabel("time (days)"),
        Guide.ylabel("S_c"), Geom.line)
    p5 = plot(solDataFrame, x=:x1, y=:x6,
        Guide.xlabel("time (days)"),
        Guide.ylabel("L_c"), Geom.line)
    p6 = plot(solDataFrame, x=:x1, y=:x7,
        Guide.xlabel("time (days)"),
        Guide.ylabel("I_cl"), Geom.line)
    p7 = plot(solDataFrame, x=:x1, y=:x8,
            Guide.xlabel("time (days)"),
            Guide.ylabel("I_ch"), Geom.line)
    p8 = plot(solDataFrame, x=:x1, y=:x9,
                Guide.xlabel("time (days)"),
                Guide.ylabel("T_c"), Geom.line)
    #
    n_flyes = solDataFrame[:, 2 : 4]
    n_cows = solDataFrame[:, 5 : 9]
    #
    time = solDataFrame[:, 1]
    n_flyes = n_flyes.x2[:] + n_flyes.x3[:] + n_flyes.x4[:]
    n_cows = n_cows.x5[:] + n_cows.x6[:] + n_cows.x7[:] +
            + n_cows.x8[:] + n_cows.x9[:]
        #
    p9 = plot(x=time, y=n_flyes,
        Guide.xlabel("time (days)"),
        Guide.ylabel("Flyes CL"), Geom.line)
    p10 = plot(x=time, y=n_cows,
        Guide.xlabel("time (days)"),
        Guide.ylabel("Cows CL"), Geom.line)
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
    img0 = SVG("flyes_disease_dynamics_rkf.svg", 19cm, 11.74289cm)
    img1 = SVG("cows_disease_dynamics_rkf.svg", 19cm, 11.74289cm)
    draw(img0, plt0)
    draw(img1, plt1)
end
