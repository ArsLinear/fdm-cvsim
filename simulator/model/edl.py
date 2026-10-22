import numpy as np

def edl_current(model, dt):
    current = model.C_dl * np.gradient(model.waveform, dt)
    return current