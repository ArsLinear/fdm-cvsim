import numpy as np

E = 1

F_CONST = 96485 # C/mol
R_CONST = 8.314 # J/(K mol)
TEMP = 298.15 # K
N = 1

ALPHA = 0.5
E_0 = 0
D_OX = 1e-5
D_RED = 1e-5

D = np.array(
    [D_OX, 0],
    [0, D_RED]
)

S = np.array(
    [1],
    [-1]
)

A = np.array(
    [np.exp((-ALPHA * N * F_CONST) / (R_CONST * TEMP) * (E - E_0))],
    [np.exp(((1 - ALPHA) * N * F_CONST) / (R_CONST * TEMP) * (E - E_0))]
)