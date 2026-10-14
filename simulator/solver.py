from scipy.sparse.linalg import spsolve

def solver(big_matrix, rhs):
    return spsolve(big_matrix, rhs)