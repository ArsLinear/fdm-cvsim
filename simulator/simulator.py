import numpy as np
from .model.edl import edl_current
from .solver import solver
from .assembly import (
    concentration_matrix,
    coefficient_matrix,
    boundary_matrix,
    big_matrix_init,
    big_matrix_refresh,
    rhs_assemble,
)

def faraday_sim(model, grid, dt):
    time_array = np.arange(len(model.waveform)) * dt
    faraday_current = np.zeros(len(model.waveform))
    
    init_boundary = np.array([model.c_ox, model.c_red])
    coefficient, surface_coefficient = coefficient_matrix(model.diffusion_matrix(), grid, dt)

    concentration = concentration_matrix(init_boundary, grid)

    big_matrix = big_matrix_init(coefficient, surface_coefficient)
    big_matrix = big_matrix_refresh(
        big_matrix, 
        surface_coefficient, 
        boundary_matrix(
            model.interfacial_reaction_rate_coefficient_matrix(model.waveform[0]), 
            grid[1]-grid[0], 
            dt
        )
    )
    big_matrix = big_matrix.tocsc()

    for t in range(1, len(model.waveform)):
        flux = boundary_matrix(
            model.interfacial_reaction_rate_coefficient_matrix(model.waveform[t]), 
            grid[1]-grid[0], 
            dt
        )
        big_matrix = big_matrix_refresh(big_matrix, surface_coefficient, flux)
        rhs = rhs_assemble(concentration, init_boundary)
        
        solution = solver(big_matrix, rhs)
        concentration = solution.reshape(len(grid), 2)
        faraday_current[t] = model.butler_volmer(model.waveform[t], concentration[0, :])

    return time_array, faraday_current

def edl_sim(model, dt):
    return edl_current(model, dt)

def current_sim(model, grid, dt):
    time_array, faraday_current = faraday_sim(model, grid, dt)
    edl_current = edl_sim(model, dt)
    return time_array, faraday_current + edl_current