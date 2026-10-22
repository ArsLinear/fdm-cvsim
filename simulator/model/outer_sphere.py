import numpy as np
from ..mesh import exponential_grid, uniform_grid

R_CONST = 8.314  # J/(mol*K)
F_CONST = 96485  # C/mol
TEMP = 298.15  # K，25 °C

class OuterSphere:
    def __init__(self, waveform, k0, alpha, E0, n, D_ox, D_red, c_ox, c_red, C_dl):
        self.waveform = waveform
        self.k0 = k0
        self.alpha = alpha
        self.E0 = E0
        self.n = n
        self.D_ox = D_ox
        self.D_red = D_red
        self.c_ox = c_ox
        self.c_red = c_red
        self.C_dl = C_dl

    def diffusion_matrix(self):
        return np.array([[self.D_ox, 0.0], [0.0, self.D_red]])

    def interfacial_reaction_rate_coefficient_matrix(self, t, concentration_surface):
        potential = self.waveform[t]
        bv_matrix = np.array([np.exp(-self.alpha * self.n * F_CONST * (potential - self.E0) / (R_CONST * TEMP)), 
                              -np.exp((1 - self.alpha) * self.n * F_CONST * (potential - self.E0) / (R_CONST * TEMP)),])
        return -self.n * F_CONST * self.k0 * (bv_matrix @ concentration_surface)
