from threading import local

import numpy as np
from scipy.linalg import solve_banded
from scipy.linalg.lapack import dgttrf, dgttrs


_cache = local()


def _backsolve(factor, rhs):
    result, _ = dgttrs(*factor, rhs.copy().reshape(-1, 1), overwrite_b=True)
    return result[:, 0]


class _DiffusionSolver:
    def __init__(self, big_matrix):
        self.matrix = big_matrix
        self.fixed = big_matrix.copy()
        self.coupling = big_matrix[0, 2:4].copy()
        self.factors = []
        self.response = []

        # 节点 1..N 的扩散部分固定，两种物质分别分解一次。
        for species in range(2):
            columns = np.arange(2 + species, big_matrix.shape[1], 2)
            *factor, info = dgttrf(
                big_matrix[4, columns[:-1]], big_matrix[2, columns], big_matrix[0, columns[1:]],
            )
            if info:
                self.factors = None
                return
            self.factors.append(factor)

            rhs = np.zeros(len(columns))
            rhs[0] = big_matrix[4, species]
            self.response.append(_backsolve(factor, rhs))

        self.correction = self.coupling * np.array([x[0] for x in self.response])

    def matches(self, big_matrix):
        return (
            self.matrix is big_matrix
            and np.array_equal(self.fixed[:, 2:], big_matrix[:, 2:])
            and np.array_equal(self.fixed[4, :2], big_matrix[4, :2])
            and big_matrix[3, 1] == 0
        )

    def solve(self, big_matrix, rhs):
        if self.factors is None:
            return solve_banded((2, 2), big_matrix, rhs, overwrite_ab=False)
        interior = [
            _backsolve(factor, rhs[2 + species::2])
            for species, factor in enumerate(self.factors)
        ]

        # 消去内部浓度后，只需求解表面的 2×2 反应块。
        surface = np.array([
            [big_matrix[2, 0] - self.correction[0], big_matrix[1, 1]],
            [big_matrix[3, 0], big_matrix[2, 1] - self.correction[1]],
        ])
        surface_concentration = np.linalg.solve(
            surface, rhs[:2] - self.coupling * np.array([x[0] for x in interior]),
        )

        solution = np.empty_like(rhs)
        solution[:2] = surface_concentration
        for species in range(2):
            solution[2 + species::2] = (
                interior[species] - self.response[species] * surface_concentration[species]
            )
        return solution


def solver(big_matrix, rhs):
    big_matrix = np.asarray_chkfinite(big_matrix)
    rhs = np.asarray_chkfinite(rhs)
    if (
        big_matrix.ndim != 2 or big_matrix.shape[0] != 5
        or big_matrix.shape[1] < 8 or big_matrix.shape[1] % 2
        or big_matrix.dtype != np.float64 or rhs.dtype != np.float64
        or rhs.shape != (big_matrix.shape[1],)
    ):
        return solve_banded((2, 2), big_matrix, rhs, overwrite_ab=False)

    cached = getattr(_cache, "solver", None)
    if cached is None or not cached.matches(big_matrix):
        if np.any(big_matrix[1, 2:]) or np.any(big_matrix[3, 1:]):
            return solve_banded((2, 2), big_matrix, rhs, overwrite_ab=False)
        cached = _DiffusionSolver(big_matrix)
        _cache.solver = cached

    return cached.solve(big_matrix, rhs)
