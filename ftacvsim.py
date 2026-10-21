import time
import numpy as np
import simulator as sim

SPACE = 0.02  # cm，数值模拟的扩散区域长度
SPACE_NUM = 80
SAMPLE_SPACE = SPACE / SPACE_NUM


start_time = time.perf_counter()

experimental_data = np.loadtxt(
    r"C:\Users\bkzeng\Code\ftacv-params\data\K3FeCN6_100mVs_a.csv",
    delimiter=",", skiprows=4, usecols=(0, 1, 2),
)
SAMPLE_RATE = 1 / np.mean(np.diff(experimental_data[:, 0]))
waveform_array = experimental_data[:, 1] - experimental_data[:, 2] * 11.54  # iR 校正，R = 11.54 Ω

model = sim.OuterSphere(
    waveform=waveform_array,
    k0=0.012643958786395187,        # cm/s；100 mV/s_a MCMC
    alpha=0.4875426181152853,
    E0=0.30327993797592623,         # V
    c_ox=1.21e-6,            # mol/cm^3
    c_red=0.0,
    D_ox=7.74e-6,            # cm^2/s
    D_red=6.61e-6,           # cm^2/s
    C_dl=4.1108204049071266e-05,     # F/cm^2
    n=1,
    dx=SAMPLE_SPACE,
    dt=1/SAMPLE_RATE
)

time_array, current_faraday = sim.faraday_sim(model)

current_non_faraday = sim.edl_sim(model)

current = current_faraday + current_non_faraday

elapsed = time.perf_counter() - start_time

print(f"模拟耗时：{elapsed:.3f} 秒")
