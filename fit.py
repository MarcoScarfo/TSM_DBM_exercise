import matplotlib.pyplot as plt
import numpy as np


def fit_isotropic_hardening(sample, model):
    # SAMPLE CURVE EXTRACTION
    sigma = sample['stress'][sample['id_Rp02']:sample['id_Rm']]
    epsilon = sample['strain'][sample['id_Rp02']:sample['id_Rm']]

    # PARAMETERS INITIALIZATION
    A_init = sample['stress'][sample['id_Rp02']] * (1 + sample['strain'][sample['id_Rp02']])/2
    B_init = sample['stress'][sample['id_Rm']]
    RMSE_min = 1e10
    model['A'] = 0
    model['B'] = 0
    model['n'] = 0
    step = 0.1

    # GRID-SEARCH ON 3 LEVELS
    for A in np.arange(0.9*A_init, 1.1*A_init, step*A_init):
        for B in np.arange(0.6*B_init, 0.8*B_init, step*B_init):
            for n in np.arange(0, 0.2, step/10):
                sigma_fitted = A + B*epsilon**n
                RMSE = np.sqrt(sum((sigma - sigma_fitted) ** 2)) / len(sigma)
                if RMSE < RMSE_min:
                    RMSE_min = RMSE
                    model['A'] = np.round(A, 3)
                    model['B'] = np.round(B, 3)
                    model['n'] = np.round(n, 4)
                    print("RMSE: {}, A: {}, B: {}, n: {}".format(np.round(RMSE, 3), np.round(A, 3), np.round(B, 3), np.round(n, 4)))

    # GET OPTIMAL CASE
    sigma_fitted_optimal = model['A'] + model['B']*epsilon**model['n']

    # PLOT OPTIMAL CASE
    plt.Figure()
    plt.plot(epsilon, sigma)
    plt.plot(epsilon, sigma_fitted_optimal)
    plt.title("isotropic hardening fitting")
    plt.show()
    return model


def fit_strain_rate_hardening(sample, model, measure):
    # SAMPLE CURVE EXTRACTION
    sigma = sample['stress'][sample['id_Rp02']:sample['id_Rm']]
    epsilon = sample['strain'][sample['id_Rp02']:sample['id_Rm']]
    d_epsilon = measure[0]

    # PARAMETERS INITIALIZATION
    RMSE_min = 1e10
    model['C'] = 0
    d_epsilon_0 = 1  # CAN ALSO BE 1
    isotropic_hardening = model['A'] + model['B']*epsilon**model['n']
    step = 0.0001

    # GRID-SEARCH ON 1 LEVEL
    for C in np.arange(-0.5, 0.5, step):
        sigma_fitted = isotropic_hardening*(1 + C*np.log(d_epsilon/d_epsilon_0))
        RMSE = np.sqrt(sum((sigma - sigma_fitted)**2))/len(sigma)
        if RMSE < RMSE_min:
            RMSE_min = RMSE
            model['C'] = np.round(C, 3)
        print("RMSE: {}, C: {}".format(np.round(RMSE, 3), np.round(C, 3)))

    # GET OPTIMAL CASE
    sigma_fitted_optimal = isotropic_hardening * (1 + model['C']*np.log(d_epsilon/d_epsilon_0))

    # PLOT OPTIMAL CASE
    plt.Figure()
    plt.plot(epsilon, sigma)
    plt.plot(epsilon, sigma_fitted_optimal)
    plt.title("strain-rate hardening fitting")
    plt.show()
    return model


def fit_thermal_softening(sample, model, measure):
    # SAMPLE CURVE EXTRACTION
    sigma = sample['stress'][sample['id_Rp02']:sample['id_Rm']]
    epsilon = sample['strain'][sample['id_Rp02']:sample['id_Rm']]
    d_epsilon = measure[0]

    # TEMPERATURES
    T = measure[1]
    T_m = 1400  # 1370-1430 °C
    T_r = 20
    if T < T_r:
        T_hat = 0
    elif (T > T_r) and (T < T_m):
        T_hat = (T - T_r) / (T_m - T_r)
    else:
        T_hat = 1

    # PARAMETERS INITIALIZATION
    RMSE_min = 1e10
    model['m'] = 0
    d_epsilon_0 = 0.001  # CAN ALSO BE 1
    isotropic_hardening = model['A'] + model['B']*epsilon**model['n']
    strain_rate_hardening = 1 + model['C']*np.log(d_epsilon/d_epsilon_0)
    step = 0.001

    # GRID-SEARCH ON 1 LEVEL
    for m in np.arange(0.5, 1, step):
        sigma_fitted = isotropic_hardening * strain_rate_hardening * (1 - T_hat**m)
        RMSE = np.sqrt(sum((sigma - sigma_fitted) ** 2)) / len(sigma)
        if RMSE < RMSE_min:
            RMSE_min = RMSE
            model['m'] = np.round(m, 4)
            print("RMSE: {}, m: {}".format(np.round(RMSE, 3), np.round(m, 4)))

    # GET OPTIMAL CASE
    sigma_fitted_optimal = isotropic_hardening * strain_rate_hardening * (1 - T_hat**model['m'])

    # PLOT OPTIMAL CASE
    plt.Figure()
    plt.plot(epsilon, sigma)
    plt.plot(epsilon, sigma_fitted_optimal)
    plt.title("Thermal softening fitting")
    plt.show()
    return model
