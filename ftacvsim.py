import time
import numpy as np
import simulator as sim

SPACE = 0.02  # cm，数值模拟的扩散区域长度
SPACE_NUM = 80
SAMPLE_SPACE = SPACE / SPACE_NUM

# Ru(NH3)6^3+/2+，玻碳电极，1 M KNO3，Ag/AgCl (3 M KCl) 参比。
# 实验条件和 k0、E0：doi.org/10.1038/s41598-024-67840-x
# D_ox、D_red：doi.org/10.1016/j.jelechem.2010.12.011
# 扩散系数测于 0.5 M KNO3；此处沿用前一篇论文将其用于 1 M KNO3 的近似。
# 单位：长度 cm，时间 s，浓度 mol/cm^3，电位 V，输出电流密度 A/cm^2。
# 时间步、扩散区域长度和网格数是数值设置，不是实验测量值。

start_time = time.perf_counter()

experimental_data = np.loadtxt(
    r"C:\Users\bkzeng\Code\ftacv-params\data\K3FeCN6_100mVs_a.csv",
    delimiter=",", skiprows=4, usecols=(0, 1, 2),
)
SAMPLE_RATE = 1 / np.mean(np.diff(experimental_data[:, 0]))
waveform_array = experimental_data[:, 1] - experimental_data[:, 2] * 11.54  # iR 校正，R = 11.54 Ω

model = sim.OuterSphere(
    waveform=waveform_array,
    k0=0.0126434249,        # cm/s；100 mV/s_a MCMC
    alpha=0.487558372,
    E0=0.303279851,         # V
    c_ox=1.21e-6,            # mol/cm^3
    c_red=0.0,
    D_ox=7.74e-6,            # cm^2/s
    D_red=6.61e-6,           # cm^2/s
    C_dl=4.11089089e-05,     # F/cm^2
    n=1,
    dx=SAMPLE_SPACE,
    dt=1/SAMPLE_RATE
)

time_array, current_faraday = sim.faraday_sim(model)

current_non_faraday = sim.edl_sim(model)

current = current_faraday + current_non_faraday

elapsed = time.perf_counter() - start_time

print(f"模拟耗时：{elapsed:.3f} 秒")
