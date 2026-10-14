from .edl import edl_current
from .outer_sphere import faraday_sim


def current_sim(waveform, k0, alpha, E0, c_ox, c_red, D_ox, D_red, n, C_dl, dx, dt):
    
    time_array, faraday_current = faraday_sim(waveform, k0, alpha, E0, c_ox, c_red, D_ox, D_red, n, dx, dt)
    
    non_faraday_current = edl_current(waveform, C_dl, dt)
    
    total_current = faraday_current + non_faraday_current
    
    return time_array, total_current