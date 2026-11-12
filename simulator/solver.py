from scipy.linalg import solve_banded


def solver(big_matrix, rhs):
    return solve_banded((2, 2), big_matrix, rhs, overwrite_ab=False)