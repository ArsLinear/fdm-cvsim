import numpy as np

SAMPLE_TIME = 0.001

SPACE = 1
SPACE_NUM = 1000
SAMPLE_SPACE = SPACE / SPACE_NUM

R_CONST = 8.314  # J/(mol*K)
F_CONST = 96485  # C/mol
TEMP = 298.15  # K

def waveform_generator(initial, switch, scan_rate):

    dt = SAMPLE_TIME

    t_1 = np.arange(0, (switch - initial) / scan_rate, dt)
    t_2 = np.arange((switch - initial) / scan_rate, 2 * (switch - initial) / scan_rate, dt)

    waveform_1 = initial + scan_rate * t_1
    waveform_2 = switch - scan_rate * (t_2 - (switch - initial) / scan_rate)

    waveform = np.concatenate([waveform_1, waveform_2])

    return waveform


def _fdm_sim(waveform, k0, alpha, E0, c_ox, c_red, D_ox, D_red, n):

    dt = SAMPLE_TIME
    dx = SAMPLE_SPACE
    num = SPACE_NUM + 1

    concentration = np.zeros((num, 2))

    for i in range(0, num):
        concentration[i, 0] = c_ox
        concentration[i, 1] = c_red

    diffusivity = np.array([[D_ox, 0.0], [0.0, D_red]])
    coefficient = diffusivity * dt / dx ** 2

    butler_volmer = np.array([np.exp(-alpha * n * F_CONST * (waveform - E0) / (R_CONST * TEMP)),
                             - np.exp((1 - alpha) * n * F_CONST * (waveform - E0) / (R_CONST * TEMP))])
    flux = 2 * k0 * dt / dx * butler_volmer

    big_martrix = np.zeros((num, num))

    big_martrix[0, 0] = np.ones(2) + coefficient - flux
    big_martrix[0, 1] = -2 * coefficient

    for i in range(1, num-1):
        big_martrix[i, i-1] = -coefficient
        big_martrix[i, i] = np.ones(2) + 2 * coefficient
        big_martrix[i, i+1] = -coefficient

    big_martrix[num, num] = np.array([c_ox, 0], [0, c_red])

    concentration = np.linalg.solve(big_martrix, concentration)

    return concentration

def cv_sim(initial, switch, scan_rate, k0, alpha, E0, c_ox, c_red, D_ox, D_red, n):

    step = int((switch - initial) / scan_rate / SAMPLE_TIME)
    current = np.zeros(step)

    waveform = waveform_generator(initial, switch, scan_rate)

    for t in range(len(waveform)):
        concentration = _fdm_sim(waveform[t], k0, alpha, E0, c_ox, c_red, D_ox, D_red, n)
        concentration_surface = concentration[0, :]
        current[t] = n * F_CONST * D_red * (concentration[1, 1] - concentration[0, 1]) / SAMPLE_SPACE
        butler_volmer = np.array([np.exp(-alpha * n * F_CONST * (waveform - E0) / (R_CONST * TEMP)),
                            - np.exp((1 - alpha) * n * F_CONST * (waveform - E0) / (R_CONST * TEMP))])
        current[t] = n * F_CONST * k0 * butler_volmer * concentration_surface

    return [waveform, current]