import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

def solver(coefficient, flux, concentration, far_boundary):
    num = concentration.shape[0]

    big_matrix = np.zeros((num, num, 2, 2))
    
    for i in range(1, num - 1):
        big_matrix[i, i-1] = -coefficient
        big_matrix[i, i] = np.eye(2) + 2 * coefficient
        big_matrix[i, i+1] = -coefficient

    big_matrix[num-1, num-1] = np.eye(2)
    big_matrix[0, 1] = -2 * coefficient
    big_matrix[0, 0] = np.eye(2) + 2 * coefficient + flux

    matrix_2d = big_matrix.transpose(0, 2, 1, 3).reshape(2 * num, 2 * num)
    concentration[num-1] = far_boundary
    concentration_matrix = spsolve(csr_matrix(matrix_2d), concentration.reshape(-1)).reshape(num, 2)

    return concentration_matrix