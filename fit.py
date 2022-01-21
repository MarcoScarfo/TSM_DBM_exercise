import matplotlib.pyplot as plt


def fit_isotropic_hardening(eps, sig):
    # PARAMETERS INITIALIZATION
    MSE_min = 1e10
    A_opt = 0
    B_opt = 0
    n_opt = 0

    # GRID-SEARCH ON 3 LEVELS
    for A in range(0, 2, 0.1):
        for B in range(0, 2, 0.1):
            for n in range(0, 2, 0.1):
                sig_fit = A + B*(eps)**n
                MSE = sum((sig - sig_fit)**2)
                if MSE < MSE_min:
                    MSE_min = MSE
                    A_opt = A;
                    B_opt = B;
                    n_opt = n;
                    print("A: {}, B: {}, n: {}".format(A, B, n))

    # GET OPTIMAL CASE
    sig_fit_opt = A_opt + B_opt*(eps)**n_opt

    # PLOT OPTIMAL CASE
    plt.Figure()
    plt.scatter(eps, sig)
    plt.plot(eps, sig_fit_opt)
    plt.show()
    return [A_opt, B_opt, n_opt]


def fit_strain_rate_hardening(eps, sig, A, B, n):
    # PARAMETERS INITIALIZATION
    MSE_min = 1e10
    C_opt = 0
    eps_0 =

    # GRID-SEARCH ON 3 LEVELS
    for C in range(0, 2, 0.1):
                sig_fit = (A + B*(eps)**n)*(1 + C*np.log(eps/eps_0))
                MSE = sum((sig - sig_fit)**2)
                if MSE < MSE_min:
                    MSE_min = MSE
                    A_opt = A;
                    B_opt = B;
                    N_opt = N;
                    print("A: {}, B: {}, N: {}".format(A, B, N))

    # GET OPTIMAL CASE
    sig_fit_opt = A_opt + B_opt*(eps)**N_opt

    # PLOT OPTIMAL CASE
    plt.Figure()
    plt.scatter(eps, sig)
    plt.plot(eps, sig_fit_opt)
    plt.show()
    return C


def fit_thermal_softening(eps, sig):
    pass