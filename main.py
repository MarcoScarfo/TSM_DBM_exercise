import pandas as pd
from processing import *
from plot import *
from fit import *
import os

# SYSTEM PARAMETERS
plot = False
model = dict()
cases = list()

# BAR PARAMETERS
bar = dict()
bar['d'] = 10e-3  # Diametro barra
bar['A'] = (bar['d'] / 2) ** 2 * np.pi  # Area barra
bar['E'] = 1.928e11  # Modulo elastico in Pa
bar['rho'] = 7.9e3  # kg/m3

# SAMPLE PARAMETERS
sample = dict()
sample['d'] = 3e-3  # Diametro provino
sample['L'] = 5e-2  # Lunghezza provino
sample['A'] = (sample['d'] / 2) ** 2 * np.pi  # Area provino
sample['E'] = 2e11  # Modulo elastico in Pa
sample['rho'] = 8.17e3  # kg/m3
sample['Co'] = np.sqrt(sample['E'] / sample['rho'])  # Velocità del suono nel provino in m/s

# RAW DATA
file = list()
file.append(r'\Data\Inconel 718\Inc718-Cast_050 - Statica - 20C\Inc718-Cast_50.txt')
file.append(r'\Data\Inconel 718\Inc718-Cast_001 - 200 s-1 - 20C\Cast-200-20C-001-in-out-signals.txt')
file.append(r'\Data\Inconel 718\Inc718-Cast_007 - 800 s-1 - 20C\Cast-800-20C-007-in-out-signals.txt')
file.append(r'\Data\Inconel 718\Inc718-Cast_017 - 200 s-1 - 550C\Cast-200-550C-017-all_data.txt')

# MEASURES DATA
measures = list()
measures.append([0, 20])  # Case Inc718-Cast_50
measures.append([200, 20])  # Case Cast-200-20C-001-in-out-signals
measures.append([800, 20])  # Case Cast-800-20C-007-in-out-signals
measures.append([200, 550])  # Case Cast-200-550C-017-all_data

# FRACTURE DATA
fracture_data = list()
fracture_data.append([2.220e-3, 6.490e-3])  # Case Inc718-Cast_50
fracture_data.append([2.416e-3, 5.155e-3])  # Case Cast-200-20C-001-in-out-signals
fracture_data.append([2.292e-3, 5.155e-3])  # Case Cast-800-20C-007-in-out-signals

# START AND END IMPULSE IDS
pulse_ids = list()
pulse_ids.append([0, 0])  # Case Cast-200-20C-001-in-out-signals
pulse_ids.append([37, 219])  # Case Cast-200-20C-001-in-out-signals
pulse_ids.append([315, 861])  # Case Cast-800-20C-007-in-out-signals
pulse_ids.append([0, 0])  # Case Cast-200-550C-017-all_data

# START AND END IMPULSE IDS
linear_ids = list()
linear_ids.append([0, 25000])  # Case Cast-200-20C-001-in-out-signals
linear_ids.append([3, 5])  # Case Cast-200-20C-001-in-out-signals
linear_ids.append([30, 40])  # Case Cast-800-20C-007-in-out-signals
linear_ids.append([0, 38])  # Case Cast-200-550C-017-all_data

# LOOP ON 4 CASES
for id_case in range(0, 4):
    # DATA READING
    file_path = os.path.dirname(os.path.realpath(__file__)) + file[id_case]
    data_frame = pd.read_csv(file_path, sep="\t", index_col=False)
    column_names = list(data_frame.columns)

    # TIME EXTRACTION
    time = data_frame['Time [s]'].values

    if (id_case == 0) or (id_case == 3):
        # STRESS AND STRAIN ALREADY DEFINED
        sample['stress'] = data_frame[column_names[2]]
        sample['strain'] = data_frame[column_names[1]]
    else:
        # FORCE TRANSMITTED IN BARS
        force = dict()
        force['transmitted'] = data_frame[column_names[2]]
        force['reflected'] = data_frame[column_names[1]]

        # STRESS, STRAIN AND STRAIN-RATE EXTRACTION
        sample, bar = from_forces_to_strain_stress(time, bar, sample, force, id_case, pulse_ids)

    # YOUNG'S MODULUS CORRECTION
    real = dict()
    # sample['diff stress'] = np.diff(sample['stress'])
    # ids = np.where(sample['diff stress'] > (np.max(sample['diff stress']) - np.min(sample['diff stress']))/2)[0]
    # id = np.argmax(sample['diff stress'])
    real['sigma'] = sample['stress'][linear_ids[id_case][1]] - sample['stress'][linear_ids[id_case][0]]
    real['epsilon'] = sample['strain'][linear_ids[id_case][1]] - sample['strain'][linear_ids[id_case][0]]
    real['E'] = real['sigma']/real['epsilon'] * 1e6

    # Correction only on linear part
    sample['strain corrected'] = list()
    for i in range(sample['strain'].shape[0]):
        if sample['strain'][i] <= real['epsilon']:
            sample['strain corrected'].append(sample['strain'][i] - sample['stress'][i] * (1/real['E'] - 1/sample['E']) * 1e6)
        else:
            sample['strain corrected'].append(sample['strain'][i] - real['sigma'] * (1 / real['E'] - 1 / sample['E']) * 1e6)

    # Correction entire curve
    sample['stress corrected'] = sample['stress']
    sample['strain corrected'] = np.array(sample['strain corrected'])

    # Rp0.2 IDENTIFICATION
    theoretical = dict()
    theoretical['sigma'] = sample['stress corrected'][linear_ids[id_case][1]] - sample['stress corrected'][linear_ids[id_case][0]]
    theoretical['epsilon'] = sample['strain corrected'][linear_ids[id_case][1]] - sample['strain corrected'][linear_ids[id_case][0]]
    theoretical['E'] = theoretical['sigma']/theoretical['epsilon']
    theoretical['stress reconstruct'] = theoretical['E'] * (sample['strain corrected'] - 0.002)
    sample['id_Rp02'] = np.where(theoretical['stress reconstruct'] > sample['stress corrected'])[0][0]

    # NECKING IDENTIFICATION
    sample['id_Rm'] = np.argmax(sample['stress corrected'])
    sample['stress no necking'] = sample['stress corrected'][:sample['id_Rm']]
    sample['strain no necking'] = sample['strain corrected'][:sample['id_Rm']]

    # TRUE CURVES CORRECTION
    sample['strain true'] = np.log(1 + sample['strain corrected'])
    sample['stress true'] = (1 + sample['strain corrected']) * sample['stress corrected']
    sample['stress true no necking'] = sample['stress true'][:sample['id_Rm']]
    sample['strain true no necking'] = sample['strain true'][:sample['id_Rm']]

    # FRACTURE VALUES
    if (id_case == 0) or (id_case == 3):
        sample['stress extended'] = sample['stress true']
        sample['strain extended'] = sample['strain true']
    else:
        a = fracture_data[id_case][0]/2
        R = fracture_data[id_case][1]

        sample['average sigma'] = bar['force reflected'][-1] / (np.pi * a**2)
        sample['stress fracture'] = sample['average sigma']/((1 + 2*R/a) * np.log(1 + a/(2*R))) / 1e6
        sample['strain fracture'] = 2 * np.log(sample['d'] / (2 * a))
        sample['stress extended'] = np.append(sample['stress true no necking'], [sample['stress fracture']])
        sample['strain extended'] = np.append(sample['strain true no necking'], [sample['strain fracture']])

    # JOHNSON-COOK MODEL FITTING
    # ISOTROPIC HARDENING
    if id_case == 0:
        model = fit_isotropic_hardening(sample, model)

    # STRAIN-RATE HARDENING
    if id_case == 1:
        model = fit_strain_rate_hardening(sample, model, measures[id_case])

    # THERMAL SOFTENING
    if id_case == 3:
        model = fit_thermal_softening(sample, model, measures[id_case])

    # PLOTS
    if plot:
        plot_data(time, bar, sample)

    store = dict()
    store['strain corrected'] = sample['strain corrected']
    store['stress corrected'] = sample['stress corrected']
    store['strain true'] = sample['strain true']
    store['stress true'] = sample['stress true']
    store['strain extended'] = sample['strain extended']
    store['stress extended'] = sample['stress extended']

    cases.append(store)

print("The Jonshon-Cook model fits the data with the following parameters:")
print("A: {}".format(model['A']))
print("B: {}".format(model['B']))
print("C: {}".format(model['C']))
print("n: {}".format(model['n']))
print("m: {}".format(model['m']))

