import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix

def solver(coefficient, flux, concentration, far_boundary):
    num = concentration.shape[0]

    big_matrix = lil_matrix((2 * num, 2 * num))
    
    for i in range(1, num - 1):
        start = 2 * i
        rows = slice(start, start + 2)

        big_matrix[rows, start-2:start] = -coefficient
        big_matrix[rows, start:start+2] = np.eye(2) + 2*coefficient
        big_matrix[rows, start+2:start+4] = -coefficient

    start = 2 * (num - 1)
    big_matrix[0:2, 0:2] = np.eye(2) + 2*coefficient + flux
    big_matrix[0:2, 2:4] = -2*coefficient
    big_matrix[start:start+2, start:start+2] = np.eye(2)

    concentration[num-1] = far_boundary
    concentration_matrix = spsolve(
        big_matrix.tocsc(), 
        concentration.reshape(-1),
    ).reshape(num, 2)

    return concentration_matrix