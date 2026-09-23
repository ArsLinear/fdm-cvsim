import numpy as np

SAMPLE_TIME = 0.005  # s；0.5 V/s 时每步扫描 2.5 mV

SPACE = 0.02  # cm，数值模拟的扩散区域长度
SPACE_NUM = 80
SAMPLE_SPACE = SPACE / SPACE_NUM

R_CONST = 8.314  # J/(mol*K)
F_CONST = 96485  # C/mol
TEMP = 293.15  # K，20 °C

def waveform_generator(initial, switch, scan_rate):

    dt = SAMPLE_TIME

    t_1 = np.arange(0, (switch - initial) / scan_rate, dt)
    t_2 = np.arange((switch - initial) / scan_rate, 2 * (switch - initial) / scan_rate, dt)

    waveform_1 = initial + scan_rate * t_1
    waveform_2 = switch - scan_rate * (t_2 - (switch - initial) / scan_rate)

    waveform = np.concatenate([waveform_1, waveform_2])

    return waveform


def _fdm_sim(waveform, k0, alpha, E0, c_ox, c_red, D_ox, D_red, n, concentration):

    dt = SAMPLE_TIME
    dx = SAMPLE_SPACE
    num = SPACE_NUM + 1

    diffusivity = np.array([[D_ox, 0.0], [0.0, D_red]])
    coefficient = diffusivity * dt / dx ** 2

    butler_volmer = np.array([np.exp(-alpha * n * F_CONST * (waveform - E0) / (R_CONST * TEMP)),
                             - np.exp((1 - alpha) * n * F_CONST * (waveform - E0) / (R_CONST * TEMP))])
    flux = 2 * k0 * dt / dx * np.outer([1.0, -1.0], butler_volmer)

    big_martrix = np.zeros((num, num, 2, 2))

    big_martrix[0, 0] = np.eye(2) + 2 * coefficient + flux
    big_martrix[0, 1] = -2 * coefficient

    for i in range(1, num - 1):
        big_martrix[i, i-1] = -coefficient
        big_martrix[i, i] = np.eye(2) + 2 * coefficient
        big_martrix[i, i+1] = -coefficient

    big_martrix[num-1, num-1] = np.eye(2)

    big_martrix = big_martrix.transpose(0, 2, 1, 3).reshape(2 * num, 2 * num)
    concentration[num-1] = [c_ox, c_red]
    concentration = np.linalg.solve(big_martrix, concentration.reshape(-1)).reshape(num, 2)

    return concentration

def cv_sim(initial, switch, scan_rate, k0, alpha, E0, c_ox, c_red, D_ox, D_red, n):

    waveform = waveform_generator(initial, switch, scan_rate)
    current = np.zeros(len(waveform))
    concentration = np.tile([c_ox, c_red], (SPACE_NUM + 1, 1))

    for t in range(len(waveform)):
        concentration = _fdm_sim(
            waveform[t], k0, alpha, E0, c_ox, c_red, D_ox, D_red, n, concentration
        )
        concentration_surface = concentration[0, :]
        butler_volmer = np.array([np.exp(-alpha * n * F_CONST * (waveform - E0) / (R_CONST * TEMP)),
                            - np.exp((1 - alpha) * n * F_CONST * (waveform - E0) / (R_CONST * TEMP))])
        current[t] = n * F_CONST * k0 * (butler_volmer[:, t] @ concentration_surface)

    return [waveform, current]

# Ru(NH3)6^3+/2+，玻碳电极，1 M KNO3，Ag/AgCl (3 M KCl) 参比。
# 实验条件和 k0、E0：doi.org/10.1038/s41598-024-67840-x
# D_ox、D_red：doi.org/10.1016/j.jelechem.2010.12.011
# 扩散系数测于 0.5 M KNO3；此处沿用前一篇论文将其用于 1 M KNO3 的近似。
# 单位：长度 cm，时间 s，浓度 mol/cm^3，电位 V，输出电流密度 A/cm^2。
# 时间步、扩散区域长度和网格数是数值设置，不是实验测量值。

waveform, current = cv_sim(
    initial=0.067,       # V vs Ag/AgCl；从 E0 + 0.2 V 开始
    switch=-0.333,      # V vs Ag/AgCl；在 E0 - 0.2 V 反向
    scan_rate=-0.05,     # V/s；文献测试的扫描速率范围包含 0.5 V/s
    k0=0.02,            # cm/s
    alpha=0.5,          # 对称电荷转移的近似
    E0=-0.133,          # V vs Ag/AgCl (3 M KCl)
    c_ox=1e-6,          # 1 mM Ru(NH3)6^3+ = 1e-6 mol/cm^3
    c_red=0.0,
    D_ox=5.3e-6,        # cm^2/s，Ru(NH3)6^3+
    D_red=7.3e-6,       # cm^2/s，Ru(NH3)6^2+
    n=1,
)

data = np.column_stack((waveform, current))
np.savetxt(
    "cv_ruhex.csv",
    data,
    delimiter=",",
    header="potential_V_vs_AgAgCl,current_density_A_cm2",
    comments="",
    fmt="%.8g",
)
print("结果已保存到 cv_ruhex.csv")
