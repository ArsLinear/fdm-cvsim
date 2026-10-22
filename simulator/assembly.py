import numpy as np
from scipy.sparse import lil_matrix
from .model.outer_sphere import R_CONST, F_CONST, TEMP

def boundary_matrix(model, t):
    butler_volmer = np.array([np.exp(-model.alpha * model.n * F_CONST * (model.waveform[t] - model.E0) / (R_CONST * TEMP)),
                    - np.exp((1 - model.alpha) * model.n * F_CONST * (model.waveform[t] - model.E0) / (R_CONST * TEMP))])
    flux = 2 * model.k0 * model.dt / model.dx * np.outer([1.0, -1.0], butler_volmer)
    return flux

def concentration_matrix(model):
    time_array = np.arange(len(model.waveform)) * model.dt
    time_tot = time_array[-1]

    length_diff = 6 * np.sqrt(max(model.D_ox, model.D_red) * time_tot)
    space_num = int(np.ceil(length_diff / model.dx))
    num = space_num + 1

    concentration = np.tile([model.c_ox, model.c_red], (space_num + 1, 1))

    return concentration

def coefficient_matrix(model):
    return np.array([[model.D_ox, 0.0], [0.0, model.D_red]]) * model.dt / model.dx ** 2

def butler_volmer(model, t, concentration_surface):
    potential = model.waveform[t]
    bv_matrix = np.array([
        np.exp(-model.alpha * model.n * F_CONST * (potential - model.E0) / (R_CONST * TEMP)),
        -np.exp((1 - model.alpha) * model.n * F_CONST * (potential - model.E0) / (R_CONST * TEMP)),
    ])
    return -model.n * F_CONST * model.k0 * (bv_matrix @ concentration_surface)



def big_matrix_init(coefficient, concentration):
    num = concentration.shape[0]

    big_matrix = lil_matrix((2 * num, 2 * num))
    
    for i in range(1, num - 1):
        start = 2 * i
        rows = slice(start, start + 2)

        big_matrix[rows, start-2:start] = -coefficient
        big_matrix[rows, start:start+2] = np.eye(2) + 2*coefficient
        big_matrix[rows, start+2:start+4] = -coefficient

    big_matrix[0:2, 0:2] = np.eye(2) + 2 * coefficient

    big_matrix[0:2, 2:4] = -2 * coefficient

    start = 2 * (num - 1)
    big_matrix[start:start + 2, start:start + 2] = np.eye(2)
    
    return big_matrix


def big_matrix_refresh(big_matrix, coefficient, flux):
    big_matrix[0:2, 0:2] = np.eye(2) + 2*coefficient + flux
    return big_matrix

def rhs_assemble(concentration, far_boundary):
    rhs = concentration.copy()
    rhs[-1] = far_boundary
    return rhs.reshape(-1)
