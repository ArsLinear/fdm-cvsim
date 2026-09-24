import time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

R_CONST = 8.314  # J/(mol*K)
F_CONST = 96485  # C/mol
TEMP = 293.15  # K，20 °C


def faraday_sim(waveform, k0, alpha, E0, c_ox, c_red, D_ox, D_red, n, dx, dt):

    start = time.perf_counter()

    time_array = np.arange(len(waveform)) * dt
    time_tot = time_array[-1]

    length_diff = 6 * np.sqrt(max(D_ox, D_red) * time_tot)
    space_num = int(np.ceil(length_diff / dx))
    num = space_num + 1

    current = np.zeros(len(waveform))
    concentration = np.tile([c_ox, c_red], (space_num + 1, 1))

    diffusivity = np.array([[D_ox, 0.0], [0.0, D_red]])
    coefficient = diffusivity * dt / dx ** 2

    big_matrix = np.zeros((num, num, 2, 2))

    for i in range(1, num - 1):
        big_matrix[i, i-1] = -coefficient
        big_matrix[i, i] = np.eye(2) + 2 * coefficient
        big_matrix[i, i+1] = -coefficient

    big_matrix[num-1, num-1] = np.eye(2)
    big_matrix[0, 1] = -2 * coefficient

    for t in range(1, len(waveform)):

        butler_volmer = np.array([np.exp(-alpha * n * F_CONST * (waveform[t] - E0) / (R_CONST * TEMP)),
                            - np.exp((1 - alpha) * n * F_CONST * (waveform[t] - E0) / (R_CONST * TEMP))])
        flux = 2 * k0 * dt / dx * np.outer([1.0, -1.0], butler_volmer)

        big_matrix[0, 0] = np.eye(2) + 2 * coefficient + flux

        matrix_2d = big_matrix.transpose(0, 2, 1, 3).reshape(2 * num, 2 * num)
        concentration[num-1] = [c_ox, c_red]
        concentration = spsolve(csr_matrix(matrix_2d), concentration.reshape(-1)).reshape(num, 2)

        concentration_surface = concentration[0, :]
        current[t] = -n * F_CONST * k0 * (butler_volmer @ concentration_surface)

    end = time.perf_counter()
    print(f"Simulation time: {end - start:.2f} s")

    return [time_array, current]

def non_faraday_sim(waveform, C_dl, dt):
    current = C_dl * np.gradient(waveform, dt)
    return current