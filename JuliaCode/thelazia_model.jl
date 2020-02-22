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
    w_f = 0.0; v_l = 0.0; v_h = 0.0;
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
    #
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
                + v_h * i_c_h
                )
    new_l_c = beta_c / n_c * s_c * i_f - (k_c + mu_c) * l_c
    #
    new_i_c_l  = ( theta * k_c * l_c
                    - beta_c_tilde / n_c * i_c_l * i_f
                    - (mu_c + v_l) * i_c_l
                )
    #
    new_i_c_h = ( (1.0 - theta ) * k_c * l_c
                    + beta_c_tilde / n_c * i_c_l * i_f
                    - mu_c * i_c_h
                    - v_h * i_c_h
    )
    new_t_c = - (rho + mu_c) * t_c + v_l * i_c_l
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

function rhs_f(x, u)
    n_c = p[1]; lambda_f = p[2]; lambda_c = p[3];
    beta_c = p[4]; beta_c_tilde = p[4];
    beta_f = p[5];  beta_f_tilde = p[6];
    k_f = p[7]; k_c = p[8]; mu_f = p[9];
    mu_c = p[10]; rho = p[11]; theta = p[12];
    #
    w_f = u[1];
    v_l = u[2];
    v_h = u[3];
    #
    s_f = x[1]; l_f = x[2]; i_f = x[3];
    s_c = x[4]; l_c = x[5];
    i_c_l = x[6]; i_c_h = x[7];
    t_c = x[8];

    new_s_f = (
                lambda_f
                    - beta_f / n_c * s_f * i_c_l
                    - beta_f_tilde / n_c * s_f * i_c_h
                    - mu_f * s_f
                    - w_f * s_f
            )
    #
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
                + v_h * i_c_h
                )
    new_l_c = beta_c / n_c * s_c * i_f - (k_c + mu_c) * l_c
    #
    new_i_c_l  = ( theta * k_c * l_c
                    - beta_c_tilde / n_c * i_c_l * i_f
                    - (mu_c + v_l) * i_c_l
                )
    #
    new_i_c_h = ( (1.0 - theta ) * k_c * l_c
                    + beta_c_tilde / n_c * i_c_l * i_f
                    - mu_c * i_c_h
                    - v_h * i_c_h
    )
    new_t_c = - (rho + mu_c) * t_c + v_l * i_c_l
    #
    dx = [
            new_s_f;
            new_l_f;
            new_i_f;
            new_s_c;
            new_l_c;
            new_i_c_l;
            new_i_c_h;
            new_t_c;
        ]
    return dx;
end

function rhs_adjoints(x, u, psi)
    n_c_inf = p[1]; lambda_f = p[2]; lambda_c = p[3];
    beta_c = p[4]; beta_c_tilde = p[4];
    beta_f = p[5];  beta_f_tilde = p[6];
    k_f = p[7]; k_c = p[8]; mu_f = p[9];
    mu_c = p[10]; rho = p[11]; theta = p[12];
    #
    w_f_t = u[1];
    v_l_t = u[2];
    v_h_t = u[3];
    #
    s_f = x[1]; l_f = x[2]; i_f = x[3];
    s_c = x[4]; l_c = x[5];
    i_c_l = x[6]; i_c_h = x[7];
    t_c = x[8];

    psi_s_f = psi[1]; psi_l_f = psi[2]; psi_i_f = psi[3];
    psi_s_c = psi[4]; psi_l_c = psi[5];
    psi_i_c_l = psi[6]; psi_i_c_h = psi[7];
    psi_t_c = psi[8];

    new_psi_s_f = (beta_f_tilde * i_c_h / n_c_inf
                    + beta_f * i_c_l / n_c_inf) * psi_l_f
                    - (mu_f + beta_f_tilde * i_c_h / n_c_inf
                    + beta_f * i_c_l / n_c_inf + w_f_t) * psi_s_f

    new_psi_l_f = k_f * psi_i_f - (k_f + mu_f + w_f_t) * psi_l_f

    new_psi_i_f = beta_c_tilde * i_c_l * psi_i_c_h / n_c_inf
                    - beta_c_tilde * i_c_l * psi_i_c_l / n_c_inf
                    - (mu_f + w_f_t) * psi_i_f + beta_c * psi_l_c*s_c / n_c_inf
                    - beta_c * psi_s_c * s_c / n_c_inf

    new_psi_s_c = beta_c * i_f * psi_l_c / n_c_inf
                    - (mu_c + beta_c * i_f / n_c_inf) * psi_s_c

    new_psi_l_c = - k_c * psi_i_c_h * (theta - 1.0) + k_c * psi_i_c_l * theta
                    - (k_c + mu_c) * psi_l_c

    new_psi_i_c_l = beta_c_tilde * i_f * psi_i_c_h / n_c_inf
                    - (mu_c + beta_c_tilde * i_f / n_c_inf + v_l_t)*psi_i_c_l
                    + beta_f * psi_l_f * s_f / n_c_inf
                    - beta_f * psi_s_f * s_f / n_c_inf

    new_psi_i_c_h = -(mu_c + v_h_t) * psi_i_c_h
                    + beta_f_tilde * psi_l_f * s_f / n_c_inf
                    - beta_f_tilde * psi_s_f * s_f / n_c_inf
                    + psi_t_c * v_h_t

    new_psi_t_c = -(mu_c + rho) * psi_t_c

    psi = [
        new_psi_s_f;
        new_psi_l_f;
        new_psi_i_f;
        new_psi_s_c;
        new_psi_l_c;
        new_psi_i_c_l;
        new_psi_i_c_h;
        new_psi_t_c;
        ]
    return psi;
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
    # Initial conditions
    s_f_zero = lambda_f / mu_f - 100.0;
    l_f_zrop = 40.0;
    i_f_zero = 60.0;
    s_c_zero = n_c - 10.0;
    l_c = 2.0;
    i_c_l = 8.0;
    i_c_h =  0.0;
    t_c = 0.0;
    # cost functional parameterJSON
    a_l_f = 1.0; a_i_f = 1.0; a_l_c = 1.0;
    a_i_c_l = 1.0; a_i_c_h = 1.0;
    b_f = 1.0; b_c_l = 1.0; b_c_h = 1.0;

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
        "theta": 0.7,
        "s_f_zero": 22400.0,
        "l_f_zero": 40.0,
        "i_f_zero": 60.0,
        "s_c_zero": 990,
        "l_c_zero": 2.0,
        "i_c_l_zero": 8.0,
        "i_c_h_zero": 0.0,
        "t_c_zero": 0.0,
        "a_l_f": 1.0,
        "a_i_f": 1.0,
        "a_l_c": 1.0,
        "a_i_c_l": 1.0,
        "a_i_c_h": 1.0,
        "b_f": 1.0,
        "b_c_l": 1.0,
        "b_c_h": 1.0
        }""");

    parameterDict = Dict(#
        "n_c" => 1000.0, "lambda_f" => 500.0, "lambda_c" => 0.9259259259259259,
        "beta_c" => 0.001, "beta_c_tilde" => 0.001, "beta_f" => 0.001,
        "beta_f_tilde" => 0.001, "k_f" => 0.07142857142857142,
        "k_c" => 0.02857142857142857, "mu_f" => 0.022222222222222223,
        "mu_c" => 0.000925925925925926, "rho" => .025, "theta" => .001,
        "s_f_zero" => 22400.0, "l_f_zero" => 40.0, "i_f_zero" => 60.0,
        "s_c_zero" => 990.0, "l_c_zero" => 2.0, "i_c_l_zero" => 8.0,
        "i_c_h_zero" => 0.0, "t_c_zero" => 0.0,
        "a_l_f" => 1.0, "a_i_f" => 1.0, "a_l_c" => 1.0,
        "a_i_c_l" => 1.0, "a_i_c_h" => 1.0,
        "b_f" => 1.0, "b_c_l" => 1.0, "b_c_h" => 1.0
        );
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
            parameterDict["theta"];
            parameterDict["s_f_zero"];
            parameterDict["l_f_zero"];
            parameterDict["i_f_zero"];
            parameterDict["s_c_zero"];
            parameterDict["l_c_zero"];
            parameterDict["i_c_l_zero"];
            parameterDict["i_c_h_zero"];
            parameterDict["t_c_zero"];
            parameterDict["a_l_f"];
            parameterDict["a_i_f"];
            parameterDict["a_l_c"];
            parameterDict["a_i_c_l"];
            parameterDict["a_i_c_h"];
            parameterDict["b_f"];
            parameterDict["b_c_l"];
            parameterDict["b_c_h"]
    ]
    return p;
end
