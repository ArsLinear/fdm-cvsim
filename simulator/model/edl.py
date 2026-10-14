import numpy as np

def edl_current(waveform, C_dl, dt):
    current = C_dl * np.gradient(waveform, dt)
    return current