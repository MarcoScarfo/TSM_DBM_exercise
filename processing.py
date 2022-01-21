import numpy as np


def from_forces_to_strain_stress(time, bar, sample, force, id_case, pulse_ids):
    # IDS FOR START AND END OF IMPULSE
    id_start = pulse_ids[id_case][0]
    id_end = pulse_ids[id_case][1]

    # IMPULSE FILTERING
    time = time[id_start:id_end]
    time = list(np.array(time) - time[0])
    bar['force transmitted'] = force['transmitted'].values[id_start:id_end]
    bar['force reflected'] = force['reflected'].values[id_start:id_end]

    # BAR STRAIN
    bar['strain transmitted'] = bar['force transmitted'] / (bar['A'] * bar['E'])  # in Pa
    bar['strain reflected'] = bar['force reflected'] / (bar['A'] * bar['E'])  # in Pa

    # SAMPLE STRESS
    sample['stress'] = (bar['A'] / sample['A']) * bar['E'] * bar['strain transmitted'] / 1e6  # in MPa

    # SAMPLE STRAIN
    dt = np.diff(time)
    integral = [0]
    for i in range(bar['strain reflected'].shape[0] - 1):
        integral.append(integral[-1] + dt[i] * bar['strain reflected'][i + 1])
    sample['strain'] = (2 * sample['Co'] / sample['L']) * np.asarray(integral)

    # SAMPLE STRAIN RATE
    sample['strain-rate'] = (2 * sample['Co'] / sample['L']) * bar['strain reflected']

    return sample, bar
