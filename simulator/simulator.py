import numpy as np
from .model.edl import edl_current
from .solver import solver

def faraday_sim(model):
    time_array = np.arange(len(model.waveform)) * model.dt
    faraday_current = np.zeros(len(model.waveform))

    concentration = model.concentration_matrix()
    coefficient = model.coefficient_matrix()
    far_boundary = np.array([model.c_ox, model.c_red])

    for t in range(1, len(model.waveform)):
        flux = model.boundary_matrix(t)
        concentration = solver(coefficient, flux, concentration, far_boundary)
        faraday_current[t] = model.butler_volmer(t, concentration[0, :])

    return time_array, faraday_current

def edl_sim(waveform, C_dl, dt):
    return edl_current(waveform, C_dl, dt)

def current_sim(model, C_dl):
    time_array, faraday_current = faraday_sim(model)
    edl_current = edl_sim(model.waveform, C_dl, model.dt)
    return time_array, faraday_current + edl_current