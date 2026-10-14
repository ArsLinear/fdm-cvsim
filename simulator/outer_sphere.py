import numpy as np

R_CONST = 8.314  # J/(mol*K)
F_CONST = 96485  # C/mol
TEMP = 293.15  # K，20 °C

class OuterSphere:
    def __init__(self, waveform, k0, alpha, E0, n, D_ox, D_red, c_ox, c_red, dx, dt):
        self.waveform = waveform
        self.k0 = k0
        self.alpha = alpha
        self.E0 = E0
        self.n = n
        self.D_ox = D_ox
        self.D_red = D_red
        self.c_ox = c_ox
        self.c_red = c_red
        self.dx = dx
        self.dt = dt

    def boundary_matrix(self, t):
        butler_volmer = np.array([np.exp(-self.alpha * self.n * F_CONST * (self.waveform[t] - self.E0) / (R_CONST * TEMP)),
                        - np.exp((1 - self.alpha) * self.n * F_CONST * (self.waveform[t] - self.E0) / (R_CONST * TEMP))])
        flux = 2 * self.k0 * self.dt / self.dx * np.outer([1.0, -1.0], butler_volmer)
        return flux

    def concentration_matrix(self):
        time_array = np.arange(len(self.waveform)) * self.dt
        time_tot = time_array[-1]

        length_diff = 6 * np.sqrt(max(self.D_ox, self.D_red) * time_tot)
        space_num = int(np.ceil(length_diff / self.dx))
        num = space_num + 1

        concentration = np.tile([self.c_ox, self.c_red], (space_num + 1, 1))

        return concentration

    def coefficient_matrix(self):
        return np.array([[self.D_ox, 0.0], [0.0, self.D_red]]) * self.dt / self.dx ** 2

    def butler_volmer(self, t, concentration_surface):
        potential = self.waveform[t]
        bv_matrix = np.array([
            np.exp(-self.alpha * self.n * F_CONST * (potential - self.E0) / (R_CONST * TEMP)),
            -np.exp((1 - self.alpha) * self.n * F_CONST * (potential - self.E0) / (R_CONST * TEMP)),
        ])
        return -self.n * F_CONST * self.k0 * (bv_matrix @ concentration_surface)

