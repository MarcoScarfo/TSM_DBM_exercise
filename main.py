import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
from fit import *
import os

# PLOT PARAMETERS
plot = True

# BAR PARAMETERS
d_bar = 10e-3  # Diametro barra
A_bar = (d_bar / 2) ** 2 * np.pi  # Area barra
E_bar = 1.928e11 # Modulo elastico in Pa
rho_bar = 7.9e3  # kg/m3

# SAMPLE PARAMETERS
d_sample = 3e-3  # Diametro provino
L_sample = 5e-2  # Lunghezza provino
A_sample = (d_sample / 2) ** 2 * np.pi  # Area provino
E_sample = 2e11  # Modulo elastico in Pa
rho_sample = 8.17e3  # kg/m3
Co_sample = np.sqrt(E_sample / rho_sample)  # Velocità del suono nel provino in m/s

# RAW DATA
file = list()
file.append(r'\Data\Inconel 718\Inc718-Cast_050 - Statica - 20C\Inc718-Cast_50.txt')
file.append(r'\Data\Inconel 718\Inc718-Cast_001 - 200 s-1 - 20C\Cast-200-20C-001-in-out-signals.txt')
file.append(r'\Data\Inconel 718\Inc718-Cast_007 - 800 s-1 - 20C\Cast-800-20C-007-in-out-signals.txt')

# SIGMA DATA (Rp0.2 AND Rm)
failure_data = list()
failure_data.append([])  # Case Inc718-Cast_50
failure_data.append([])  # Case Cast-200-20C-001-in-out-signals
failure_data.append([])  # Case Cast-800-20C-007-in-out-signals

# FAILURE DATA
failure_data = list()
failure_data.append([2.220e-3, 6.490e-3])  # Case Inc718-Cast_50
failure_data.append([2.416e-3, 5.155e-3])  # Case Cast-200-20C-001-in-out-signals
failure_data.append([2.292e-3, 5.155e-3])  # Case Cast-800-20C-007-in-out-signals

# START AND END IMPULS IDS
pulse_ids = list()
pulse_ids.append([0, 0])  # Case Cast-200-20C-001-in-out-signals
pulse_ids.append([37, 219])  # Case Cast-200-20C-001-in-out-signals
pulse_ids.append([315, 861])  # Case Cast-800-20C-007-in-out-signals

# LOOP ON 3 CASES
for id_case in range(1, 3):
    # DATA READING
    file_path = os.path.dirname(os.path.realpath(__file__)) + file[id_case]
    data_frame = pd.read_csv(file_path, sep="\t", index_col=False)
    column_names = list(data_frame.columns)

    # TIME EXTRACTION
    time = data_frame['Time [s]'].values

    # FORCE TRANSMITTED IN BARS
    data_frame['Force transmitted'] = data_frame[column_names[2]]
    data_frame['Force reflected'] = data_frame[column_names[1]]

    # IDS FOR START AND END OF IMPULSE
    id_start = pulse_ids[id_case][0]
    id_end = pulse_ids[id_case][1]

    # IMPULSE FILTERING
    time = time[id_start:id_end]
    time = list(np.array(time) - time[0])
    bar_force_transmitted = data_frame['Force transmitted'].values[id_start:id_end]
    bar_force_reflected = data_frame['Force reflected'].values[id_start:id_end]

    # BAR STRAIN
    bar_strain_transmitted = bar_force_transmitted / (A_bar * E_bar)  # in Pa
    bar_strain_reflected = bar_force_reflected / (A_bar * E_bar)  # in Pa

    # SAMPLE STRESS
    sample_stress = (A_bar / A_sample) * E_bar * bar_strain_transmitted  # in Pa

    # SAMPLE STRAIN
    dt = np.diff(time)
    integral = [0]
    for i in range(bar_strain_reflected.shape[0]-1):
        integral.append(integral[-1] + dt[i] * bar_strain_reflected[i+1])
    sample_strain = (2 * Co_sample / L_sample) * np.asarray(integral)

    # SAMPLE STRAIN RATE
    sample_strain_rate = (2 * Co_sample / L_sample) * bar_strain_reflected

    # YOUNG'S MODULUS CORRECTION
    diff_sample_stress = np.diff(sample_stress)
    ids = np.where(diff_sample_stress > (np.max(diff_sample_stress) - np.min(diff_sample_stress))/2)[0]
    id = np.argmax(diff_sample_stress)
    sigma_real = sample_stress[id]
    epsilon_real = sample_strain[id]
    E_real = sigma_real/epsilon_real

    # Correction only on linear part
    sample_strain_corrected = list()
    for i in range(sample_strain.shape[0]):
        if sample_strain[i] <= epsilon_real:
            sample_strain_corrected.append(sample_strain[i] - sample_stress[i] * (1/E_real - 1/E_sample))
        else:
            sample_strain_corrected.append(sample_strain[i] - sigma_real * (1 / E_real - 1 / E_sample))

    # Correction entire curve
    # sample_strain_corrected = sample_strain - sample_stress * (1 / E_real - 1 / E_sample)
    sample_stress_corrected = sample_stress
    sample_strain_corrected = np.array(sample_strain_corrected)
    # Rp0.2 IDENTIFICATION
    dE = np.diff(sample_stress_corrected/sample_strain_corrected)
    id_linear_end = np.int(np.where(dE < 0)[0][0])
    id_linear_start = np.int(np.round(id_linear_end*0.2))
    E_theoretical = (sample_stress_corrected[id_linear_end] - sample_stress_corrected[id_linear_start]) /\
                    (sample_strain_corrected[id_linear_end] - sample_strain_corrected[id_linear_start])
    #theoretical_stain_reconstruct = [0.002] + sample_strain_corrected[sample_strain_corrected > 0.002]
    theoretical_stress_reconstruct = E_theoretical * (sample_strain_corrected - 0.002)
    id_Rp02 = np.where(theoretical_stress_reconstruct > sample_stress_corrected)[0][0]

    # NECKING IDENTIFICATION
    id_necking = np.argmax(sample_stress_corrected)
    sample_stress_no_necking = sample_stress_corrected[:id_necking]
    sample_strain_no_necking = sample_strain_corrected[:id_necking]

    # FAILURE VALUES
    a = failure_data[id_case][0]/2
    R = failure_data[id_case][1]

    sample_sigma_avg = bar_force_reflected[-1] / (np.pi * a**2)
    sample_stress_fracture = sample_sigma_avg/((1 + 2*R/a) * np.log(1 + a/(2*R)))
    sample_strain_fracture = 2 * np.log(d_sample / (2 * a))
    sample_stress_extended = np.append(sample_stress_no_necking, [sample_stress_fracture])
    sample_strain_extended = np.append(sample_strain_no_necking, [sample_strain_fracture])

    # TRUE CURVES CORRECTION
    sample_strain_true = np.log(1 + sample_strain_extended)
    sample_stress_true = (1 + sample_strain_extended) * sample_stress_extended

    # JOHNSONS-COOK MODEL FITTING
    # ISOTROPIC HARDENING
    if id_case == 0:
        A, B, n = fit_isotropic_hardening()
    # STRAIN-RATE HARDENING
    if id_case == 1:
        C = fit_starin_rate_hardening()

    # THERMAL SOFTENING
    if id_case == 2:
        m = fit_thermal_softening()

    # PLOTS
    if plot:
        # BAR FORCE TRANSMITTED
        fig = px.line(x=time, y=bar_force_transmitted, title='Force transmitted')
        fig.show()

        # BAR FORCE REFLECTED
        fig = px.line(x=time, y=bar_force_reflected, title='Force reflected')
        fig.show()

        # SAMPLE STRESS
        fig = px.line(x=time, y=sample_stress/1e6, title='Sample Stress')
        fig.show()

        # SAMPLE STRAIN
        fig = px.line(x=time, y=sample_strain, title='Sample Strain')
        fig.show()

        # SAMPLE STRAIN RATE
        fig = px.line(x=time, y=sample_strain_rate, title='Sample Strain Rate')
        fig.show()

        # STRESS-STRAIN ENGINEERING NOT CORRECTED
        fig = px.line(x=sample_strain, y=sample_stress/1e6, title='Stress-Strain Engineering not corrected')
        fig.show()

        # STRESS-STRAIN ENGINEERING CORRECTED
        fig = px.line(x=sample_strain_extended, y=sample_stress_extended/1e6, title='Stress-Strain Engineering corrected')
        fig.show()

        # STRESS-STRAIN TRUE
        fig = px.line(x=sample_strain_true, y=sample_stress_true/1e6, title='Stress-Strain True')
        fig.show()

    a = 1
