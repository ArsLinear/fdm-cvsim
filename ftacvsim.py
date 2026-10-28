import time
import numpy as np
import simulator as sim
from simulator.mesh import uniform_grid, exponential_grid

start_time = time.perf_counter()

experimental_data = np.loadtxt(
    r"C:\Users\bkzeng\OneDrive\01_Workplace\AY202627_Sem_1\2_ProcessedData\20261013\K3FeCN6_10mVs.csv",
    delimiter=",", skiprows=4, usecols=(0, 1, 2),
)
SAMPLE_RATE = 1 / np.mean(np.diff(experimental_data[:, 0]))
waveform_array = experimental_data[:, 1] - experimental_data[:, 2] * 11.54  # iR 校正，R = 11.54 Ω

model = sim.OuterSphere(
    waveform=waveform_array,
    k0=0.01,        # cm/s；100 mV/s_a MCMC
    alpha=0.8,
    E0=0.3,         # V
    c_ox=1e-6,            # mol/cm^3
    c_red=0.0,
    D_ox=7e-6,            # cm^2/s
    D_red=7e-6,           # cm^2/s
    C_dl=4e-05,     # F/cm^2
    n=1,
)

length = 6 * np.sqrt(max(model.D_ox, model.D_red) * (len(waveform_array) - 1) * 1/SAMPLE_RATE)  # cm，扩散层长度，6倍扩散长度

time_array, current_faraday = sim.faraday_sim(
    model, 
    exponential_grid(length, 1e-5, ratio=1.01), 
    # uniform_grid(length, 1e-5),
    1/SAMPLE_RATE
)

current_non_faraday = sim.edl_sim(model, 1/SAMPLE_RATE)

current = current_faraday + current_non_faraday

elapsed = time.perf_counter() - start_time

print(f"模拟耗时：{elapsed:.3f} 秒")
