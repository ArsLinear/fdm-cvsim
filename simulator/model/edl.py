import numpy as np

def edl_current(model):
    current = model.C_dl * np.gradient(model.waveform, model.dt)
    return current