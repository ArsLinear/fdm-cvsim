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


def _fdm_sim(waveform, k0, alpha, E0, c_ox, c_red, D_ox, D_red, n, concentration):

    dt = SAMPLE_TIME
    dx = SAMPLE_SPACE
    num = SPACE_NUM + 1

    diffusivity = np.array([[D_ox, 0.0], [0.0, D_red]])
    coefficient = diffusivity * dt / dx ** 2

    butler_volmer = np.array([np.exp(-alpha * n * F_CONST * (waveform - E0) / (R_CONST * TEMP)),
                             - np.exp((1 - alpha) * n * F_CONST * (waveform - E0) / (R_CONST * TEMP))])
    flux = 2 * k0 * dt / dx * np.outer([1.0, -1.0], butler_volmer)

    big_martrix = np.zeros((num, num, 2, 2))

    big_martrix[0, 0] = np.eye(2) + 2 * coefficient + flux
    big_martrix[0, 1] = -2 * coefficient

    for i in range(1, num - 1):
        big_martrix[i, i-1] = -coefficient
        big_martrix[i, i] = np.eye(2) + 2 * coefficient
        big_martrix[i, i+1] = -coefficient

    big_martrix[num-1, num-1] = np.eye(2)

    big_martrix = big_martrix.transpose(0, 2, 1, 3).reshape(2 * num, 2 * num)
    concentration[num-1] = [c_ox, c_red]
    concentration = np.linalg.solve(big_martrix, concentration.reshape(-1)).reshape(num, 2)

    return concentration

def cv_sim(initial, switch, scan_rate, k0, alpha, E0, c_ox, c_red, D_ox, D_red, n):

    waveform = waveform_generator(initial, switch, scan_rate)
    current = np.zeros(len(waveform))
    concentration = np.tile([c_ox, c_red], (SPACE_NUM + 1, 1))

    for t in range(len(waveform)):
        concentration = _fdm_sim(
            waveform[t], k0, alpha, E0, c_ox, c_red, D_ox, D_red, n, concentration
        )
        concentration_surface = concentration[0, :]
        butler_volmer = np.array([np.exp(-alpha * n * F_CONST * (waveform - E0) / (R_CONST * TEMP)),
                            - np.exp((1 - alpha) * n * F_CONST * (waveform - E0) / (R_CONST * TEMP))])
        current[t] = n * F_CONST * k0 * (butler_volmer[:, t] @ concentration_surface)

    return [waveform, current]

waveform, current = cv_sim(
    initial=0.0,
    switch=0.2,
    scan_rate=0.05,
    k0=0.01,
    alpha=0.5,
    E0=0.1,
    c_ox=1.0,
    c_red=0.0,
    D_ox=1e-3,
    D_red=1e-3,
    n=1,
)

data = np.column_stack((waveform, current))
np.savetxt(
    "cv_output.csv",
    data,
    delimiter=",",
    header="potential,current",
    comments="",
    fmt="%.8g",
)
print("结果已保存到 cv_output.csv")
