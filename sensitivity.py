import numpy as np
import simulator as sim
from simulator.mesh import exponential_grid
import ftacvsim as base
import copy

DT = 1.0 / base.SAMPLE_RATE
LENGTH = base.length
RATIO = 1.01
PARAM_NAMES = ("k0", "alpha", "E0", "C_dl")

DX0_REF = 6.25e-7  # cm，暂定的参考网格

DX0_CANDIDATES = [
    5e-3,
    2e-3,
    1e-3,
    5e-4,
    2e-4,
    1e-4,
    5e-5,
    2e-5,
    1e-5,
    5e-6,
    2e-6,
    1e-6,
]

PARAM_TOL = 5e-3

def grid(dx0):
    grid = exponential_grid(LENGTH, dx0, ratio=RATIO)
    grid = grid[grid < LENGTH]
    grid = np.append(grid, LENGTH)
    return grid

def sim_current(dx0, model):
    _, current = sim.current_sim(model, grid(dx0), DT,)
    return current

def model_with_parameter(name, value):
    model = copy.copy(base.model)
    setattr(model, name, value)
    return model

def calculate_jacobian(dx0_ref):
    reference_current = sim_current(dx0_ref, base.model)

    jacobian = np.empty((len(reference_current), len(PARAM_NAMES)))

    for j, name in enumerate(PARAM_NAMES):
        eta = 1e-4
        
        value = getattr(base.model, name)

        difference_current_matrix = np.array([sim_current(dx0_ref, model_with_parameter(name, value * (1 + eta))),
                                    sim_current(dx0_ref, model_with_parameter(name, value * (1 - eta)))])
        jacobian[:, j] = np.dot([1, -1], difference_current_matrix) / (2 * eta)

    return reference_current, jacobian

reference_current, jacobian = calculate_jacobian(DX0_REF)

def estimate_parameter_shift(dx0, reference_current, jacobian):
    current = sim_current(dx0, base.model)
    delta_z, *_ = np.linalg.lstsq(jacobian, -(current - reference_current), rcond=None,)
    return delta_z

for dx0 in DX0_CANDIDATES:
    delta_z = estimate_parameter_shift(dx0, reference_current, jacobian)

    max_shift = np.max(np.abs(delta_z))
    status = "PASS" if max_shift <= PARAM_TOL else "FAIL"

    print(f"\ndx0 = {dx0:.3e} cm | {status}")

    for name, shift in zip(PARAM_NAMES, delta_z):
        print(f"  {name:>6s}: {shift:+.4e} ({shift * 100:+.4f}%)")
