import numpy as np
import waveform as wf
import simulator as sim
import time
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

SAMPLE_TIME = 0.005  # s；0.5 V/s 时每步扫描 2.5 mV

SPACE = 0.02  # cm，数值模拟的扩散区域长度
SPACE_NUM = 80
SAMPLE_SPACE = SPACE / SPACE_NUM

R_CONST = 8.314  # J/(mol*K)
F_CONST = 96485  # C/mol
TEMP = 293.15  # K，20 °C

# Ru(NH3)6^3+/2+，玻碳电极，1 M KNO3，Ag/AgCl (3 M KCl) 参比。
# 实验条件和 k0、E0：doi.org/10.1038/s41598-024-67840-x
# D_ox、D_red：doi.org/10.1016/j.jelechem.2010.12.011
# 扩散系数测于 0.5 M KNO3；此处沿用前一篇论文将其用于 1 M KNO3 的近似。
# 单位：长度 cm，时间 s，浓度 mol/cm^3，电位 V，输出电流密度 A/cm^2。
# 时间步、扩散区域长度和网格数是数值设置，不是实验测量值。

time_array, current_faraday = sim.faraday_sim(
    waveform=wf.ftacv_waveform_generator(
        initial=0.067,       # V vs Ag/AgCl；从 E0 + 0.2 V 开始
        switch=-0.333,      # V vs Ag/AgCl；在 E0 - 0.2 V 反向
        scan_rate=-0.05,     # V/s；文献测试的扫描速率范围包含 0.5 V/s
        amplitude=0.1,       # V；小幅度正弦波
        frequency=9.0,        # Hz；正弦波频率
        dt=SAMPLE_TIME
        ),
    k0=0.02,            # cm/s
    alpha=0.5,          # 对称电荷转移的近似
    E0=-0.133,          # V vs Ag/AgCl (3 M KCl)
    c_ox=1e-6,          # 1 mM Ru(NH3)6^3+ = 1e-6 mol/cm^3
    c_red=0.0,
    D_ox=5.3e-6,        # cm^2/s，Ru(NH3)6^3+
    D_red=7.3e-6,       # cm^2/s，Ru(NH3)6^2+
    n=1,
    dx=SAMPLE_SPACE,
    dt=SAMPLE_TIME
)

current_non_faraday = sim.non_faraday_sim(
    waveform=wf.ftacv_waveform_generator(
        initial=0.067,       # V vs Ag/AgCl；从 E0 + 0.2 V 开始
        switch=-0.333,      # V vs Ag/AgCl；在 E0 - 0.2 V 反向
        scan_rate=-0.05,     # V/s；文献测试的扫描速率范围包含 0.5 V/s
        amplitude=0.1,       # V；小幅度正弦波
        frequency=9.0,        # Hz；正弦波频率
        dt=SAMPLE_TIME
        ), 
    C_dl=1e-4,
    dt=SAMPLE_TIME
    )  # F/cm^2，双电层电容

current = current_faraday + current_non_faraday