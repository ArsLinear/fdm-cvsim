import numpy as np
from .model.edl import edl_current
from .solver import solver
from .assembly import big_matrix_init, big_matrix_refresh, rhs_assemble

def faraday_sim(model):
    time_array = np.arange(len(model.waveform)) * model.dt
    faraday_current = np.zeros(len(model.waveform))

    concentration = model.concentration_matrix()
    coefficient = model.coefficient_matrix()
    far_boundary = np.array([model.c_ox, model.c_red])

    num = concentration.shape[0]

    big_matrix = big_matrix_init(coefficient, concentration)
    big_matrix = big_matrix_refresh(
        big_matrix, coefficient, model.boundary_matrix(0)
    )
    big_matrix = big_matrix.tocsc()

    for t in range(1, len(model.waveform)):
        flux = model.boundary_matrix(t)
        big_matrix = big_matrix_refresh(big_matrix, coefficient, flux)
        rhs = rhs_assemble(concentration, far_boundary)
        
        solution = solver(big_matrix, rhs)
        concentration = solution.reshape(num, 2)
        faraday_current[t] = model.butler_volmer(t, concentration[0, :])

    return time_array, faraday_current

def edl_sim(model):
    return edl_current(model)

def current_sim(model):
    time_array, faraday_current = faraday_sim(model)
    edl_current = edl_sim(model)
    return time_array, faraday_current + edl_current