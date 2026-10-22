import numpy as np

def exponential_grid(length, dx0, ratio):
    num_intervals = int(np.ceil(np.log1p((ratio - 1.0) * length / dx0) / np.log1p(ratio - 1.0)))
    grid = np.concatenate(([0.0], np.cumsum(dx0 * ratio ** np.arange(num_intervals))))
    return grid

def uniform_grid(length, dx):
    num_intervals = int(np.ceil(length / dx))
    grid = np.arange(num_intervals + 1) * dx
    return grid