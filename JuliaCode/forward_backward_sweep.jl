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
        k_2 = rhs_f(x_i + 0.5 * h * k_1, u_i)
        k_3 = rhs_f(x_i + 0.5 * h * k_2, u_mean)
        k_4 = rhs_f(x_i + h * k_3, u_next)

        x_next = h / 6.0 * (k_1 + 2 * k_2 + 2 * k_3 + k_4)
        x[i + 1, :] = x_next
    end
    return x
end

function runge_kutta_backward(x, u)
    psi_s_f_final = 0.0;
    psi_l_f_final = 0.0;
    psi_i_f_final = 0.0;
    psi_s_c_final = 0.0;
    psi_l_c_final = 0.0;
    psi_i_c_l_final = 0.0;
    psi_i_c_h_final = 0.0;
    psi_t_c_final = 0.0;
    psi_final = [
        psi_s_f_final;
        psi_l_f_final;
        psi_i_f_final;
        psi_s_c_final;
        psi_l_c_final;
        psi_i_c_l_final;
        psi_i_c_h_final;
        psi_t_c_final;
        ];
    psi = zeros(n_max, x_dim)
    psi[n_max, :] = psi_final
    #
    for i in n_max : -1: 2
        psi_i = psi[i, :]
        u_i = u[i, :]
        u_previous = u[i - 1, :]
        u_mean = 0.5 * (u_i + u_previous)
        x_i = x[i, :]
        x_previous = x[i - 1, :]
        x_mean = 0.5 * (x_i + x_previous)
        #
        k_1 = rhs_adjoints(x_i, u_i, psi_i)
        k_2 = rhs_adjoints(x_mean, u_mean, psi_i - 0.5 * h * k_1)
        k_3 = rhs_adjoints(x_mean, u_mean, psi_i - 0.5 * h * k_2)
        k_4 = rhs_adjoints(x_previous, u_previous, psi_i - h * k_3 )
        #
        psi_previous = h / 6.0 * (k_1 + 2 * k_2 + 2 * k_3 + k_4)
        psi[i - 1, :] = psi_previous
    end
    return psi
end
function optimality_condition(x, psi)
    s_f = x[:, 1]; l_f = x[:, 2]; i_f = x[:, 3];
    s_c = x[:, 4]; l_c = x[:, 5];
    i_c_l = x[:, 6]; i_c_h = x[:, 7];
    t_c = x[:, 8];

    psi_s_f = psi[:, 1]; psi_l_f = psi[:, 2]; psi_i_f = psi[:, 3];
    psi_s_c = psi[:, 4]; psi_l_c = psi[:, 5];
    psi_i_c_l = psi[:, 6]; psi_i_c_h = psi[:, 7];
    psi_t_c = psi[:, 8];

    a_l_f = p[22]; a_i_f = p[23]; a_l_c = p[24];
    a_i_c_l = p[25]; a_i_c_h = p[26];
    b_f = p[27]; b_c_l = p[28]; b_c_h = p[29];

    w_f_t_aster = 0.5 * (i_f .* psi_i_f + l_f .* psi_l_f + psi_s_f .* s_f) / b_f
    v_l_t_aster = 0.5 * i_c_l .* psi_i_c_l / b_c_l
    v_h_t_aster = 0.5 * (i_c_h .* psi_i_c_h - i_c_h .* psi_t_c) / b_c_h
    # Bouds of control signals

    w_f_t_aster = min(
                        max(w_f_min * ones(n_max), w_f_t_aster),
                         w_f_max * ones(n_max)
                    )
    v_l_t_aster = min(
                        max(v_l_min * ones(n_max), v_l_t_aster),
                         v_l_max * ones(n_max)
                    )
    v_h_t_aster = min(
                        max(v_h_min * ones(n_max), v_h_t_aster),
                        v_h_max * ones(n_max)
                    )
    u_new = zeros((n_max, u_dim))
    u_new[:, 1] = w_f_t_aster
    u_new[:, 2] = v_l_t_aster
    u_new[:, 3] = v_h_t_aster
    return u_new;
end


function forward_plot()
    u_control = zeros(n_max)
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
            hstack(p8, p10)_contro
            )
    img0 = SVG("flyes_disease_dynamics_rkf.svg", 19cm, 11.74289cm)
    img1 = SVG("cows_disease_dynamics_rkf.svg", 19cm, 11.74289cm)
    draw(img0, plt0)
    draw(img1, plt1)
end

function backward_plot()
    u_control = zeros(n_max)
    x_path = runge_kutta_forward(u_control)
    psi_path = runge_kutta_backward(x_path, u_control)
    psiDataFrame=[t_span, psi_path[:, 1], psi_path[:, 2], psi_path[:, 3],
                    psi_path[:, 4], psi_path[:, 5], psi_path[:, 6],
                    psi_path[:, 7], psi_path[:, 8]]
    psiDataFrame = DataFrame(psiDataFrame)
    p1 = plot(psiDataFrame, x=:x1, y=:x2,
        Guide.xlabel("time (days)"),
        Guide.ylabel("S_f"), Geom.line)
    p2 = plot(psiDataFrame, x=:x1, y=:x3,
        Guide.xlabel("time (days)"),
        Guide.ylabel("L_f"), Geom.line)
    p3 = plot(psiDataFrame, x=:x1, y=:x4,
        Guide.xlabel("time (days)"),
        Guide.ylabel("I_f"), Geom.line)
    #
    p4 = plot(psiDataFrame, x=:x1, y=:x5,
        Guide.xlabel("time (days)"),
        Guide.ylabel("S_c"), Geom.line)
    p5 = plot(psiDataFrame, x=:x1, y=:x6,
        Guide.xlabel("time (days)"),
        Guide.ylabel("L_c"), Geom.line)
    p6 = plot(psiDataFrame, x=:x1, y=:x7,
        Guide.xlabel("time (days)"),
        Guide.ylabel("I_cl"), Geom.line)
    p7 = plot(psiDataFrame, x=:x1, y=:x8,
            Guide.xlabel("time (days)"),
            Guide.ylabel("I_ch"), Geom.line)
    p8 = plot(psiDataFrame, x=:x1, y=:x9,
                Guide.xlabel("time (days)"),
                Guide.ylabel("T_c"), Geom.line)
    #
    #
    time = psiDataFrame[:, 1]
        #
    title(hstack(p1, p2, p3), "Flyes")
    title(hstack(p4, p5, p6, p7, p8), "Cows")
    plt0 = vstack(
            hstack(p1, p2, p3),
            )
    plt1 = vstack(
            # hstack(p1, p2, p3),
            hstack(p4, p5),
            hstack(p6, p7),
            hstack(p8)
            )
    img0 = SVG("flyes_disease_dynamics_adjoints.svg", 19cm, 11.74289cm)
    img1 = SVG("cows_disease_dynamics_adjoints.svg", 19cm, 11.74289cm)
    draw(img0, plt0)
    draw(img1, plt1)
end

################################################################################
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
n_max = 10000; t_f = 2000.0;
t_span = range(0.0, t_f, length=n_max)
h = t_span[2]; eps = 1e-3;
x_dim = 8; u_dim = 3;
#
u_new = zeros((n_max, u_dim))
x_new = zeros((n_max, x_dim))
psi_new = zeros((n_max, x_dim))
#
u_old = zeros((n_max, u_dim))
x_old = zeros((n_max, x_dim))
psi_old = zeros((n_max, x_dim))
#
w_f_min = 0.0; w_f_max = 1.0
v_l_min = 0.0; v_l_max = 1.0
v_h_min = 0.0; v_h_max = 1.0
#
path = string(pwd(), "/default_parameters.json")
p = load_parameters(path);
#
eps_test = -1.0
while eps_test < 0
    u_old = u_new
    x_old = x_new
    psi_old = psi_new
    x_new = runge_kutta_forward(u_old)
    psi_new = runge_kutta_backward(x_new, u_old)
    u_new = optimality_condition(x_new, psi_new)
    u = 0.5 * (u_old + u_new)
    eps_1 = h * sum(abs(u)) - sum(abs(u_old -u));
    eps_2 = h * sum(abs(x_new)) - sum(abs(x_old - x_new));
    eps_2 = h * sum(abs(psi_new)) - sum(abs(psi_old - psi_new));
    eps_test = min(eps_1, eps_2, eps_3)
end
