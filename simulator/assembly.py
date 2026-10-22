import numpy as np
from scipy.sparse import lil_matrix
from .model.outer_sphere import R_CONST, F_CONST, TEMP

def boundary_matrix(interfacial_reaction_rate_coefficient_matrix, dx0, dt):
    return 2 * dt / dx0 * interfacial_reaction_rate_coefficient_matrix

def concentration_matrix(initial_concentration, grid):
    return np.tile(initial_concentration, (len(grid), 1))

def coefficient_matrix(diffusion_matrix, grid, dt):
    num_inner = len(grid) - 2
    coefficient_matrix = np.zeros((2, num_inner, 2, 2)) #left: j=0, right: j=1

    spacing = np.diff(grid)

    for i in range(1, len(grid) - 1):
        for j in range(2):
            coefficient_matrix[j, i - 1] = (2 * diffusion_matrix * dt / (spacing[i - (1 - j)] * (spacing[i - 1] + spacing[i])))

    surface_coefficient_matrix = diffusion_matrix * dt / spacing[0] ** 2

    return coefficient_matrix, surface_coefficient_matrix


def big_matrix_init(coefficient, surface_coefficient):
    num = coefficient.shape[1] + 2

    big_matrix = lil_matrix((2 * num, 2 * num))
    
    for i in range(1, num - 1):
        start = 2 * i
        rows = slice(start, start + 2)

        big_matrix[rows, start-2:start] = -coefficient[0, i - 1]
        big_matrix[rows, start:start+2] = np.eye(2) + coefficient[0, i - 1] + coefficient[1, i - 1]
        big_matrix[rows, start+2:start+4] = -coefficient[1, i - 1]

    big_matrix[0:2, 0:2] = np.eye(2) + 2 * surface_coefficient

    big_matrix[0:2, 2:4] = -2 * surface_coefficient

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
