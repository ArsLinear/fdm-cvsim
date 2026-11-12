# fdm-cvsim 函数与接口

> 本文档由GPT-6.1 Sol辅助整理。

## 输入、单位与调用关系

设 $E$ 为传入模型的界面电位（V），$\Delta t$ 为等时间步长（s），$x$ 为从电极表面指向溶液的距离（cm）。以 $N_t$ 和 $N_x$ 分别表示电位样点数和空间节点数。浓度用 mol/cm³，扩散系数用 cm²/s，标准速率常数用 cm/s，单位面积双电层电容用 F/cm²。模拟器返回电流密度（A/cm²）。示例 CSV 的原始电流为 A，用来计算 $iR$ 降。

示例脚本先准备电位和空间网格，再创建 `OuterSphere` 对象。它分别调用 `faraday_sim` 和 `edl_sim`，把两项电流密度逐点相加。调用者也可以把模型、网格和时间步长传给 `current_sim`，直接取得合成结果。

## 入口与波形

### `ftacvsim_example.py`

模拟步骤写在模块顶层。脚本读入实验 CSV 的前三列，依次为时间 $t_i$（s）、电位 $V_i$（V）和电流 $I_i$（A）。它用相邻时间差的平均值作为 $\Delta t$。设 $R_{\mathrm u}$ 为未补偿电阻，示例取 $10\ \Omega$，模型电位为 $E_i=V_i-I_iR_{\mathrm u}$。取两种物种中较大的扩散系数 $D_{\max}=\max(D_{\mathrm{ox}},D_{\mathrm{red}})$，空间网格至少覆盖长度 $L=6\sqrt{D_{\max}(N_t-1)\Delta t}$，第一格宽度为 $10^{-5}$ cm，后续格宽按 1.01 倍增长。[源码](https://github.com/ArsLinear/fdm-cvsim/blob/9e4f6b8755270f741ec8985a9a6233e1cde551ea/ftacvsim_example.py)

| 接口 | 输入与输出 | 实现与调用 |
| --- | --- | --- |
| 脚本顶层 | 英文逗号分隔的数值 CSV，前 1 行为表头。执行后，模块提供 `time_array`、`waveform_array`、`current_faraday`、`current_non_faraday` 和 `current`。另提供采样率 `SAMPLE_RATE`、模型 `model` 和网格目标长度 `length`。 | 将 `path/to/your/datafile.csv` 替换为实际路径。`np.loadtxt` 读取数据，`OuterSphere(...)` 保存模型参数。`exponential_grid` 建网格，`faraday_sim` 和 `edl_sim` 分别计算两项电流密度。`current` 是两项之和。顶层代码执行或导入时运行，只打印耗时，结果数组留在模块中。 |

CSV 应至少有两个样点，时间严格递增且等间隔。示例用平均时间差生成等间隔时间轴，返回的 `time_array` 从零开始。修改实验条件时，应核对模型参数、未补偿电阻和表头行数。若要使用后面的导出或灵敏度脚本，可把配置好的示例保存为 `ftacvsim.py`。仓库的 `.gitignore` 已列出这个本地文件名。[.gitignore](https://github.com/ArsLinear/fdm-cvsim/blob/9e4f6b8755270f741ec8985a9a6233e1cde551ea/.gitignore)

### `waveform.py`

设 $E_{\mathrm{initial}}$ 和 $E_{\mathrm{switch}}$ 为起始电位和转折电位，$v$ 为第一段扫描的带符号扫速。单程时间为 $t_{\mathrm{sw}}=(E_{\mathrm{switch}}-E_{\mathrm{initial}})/v$，要求 $t_{\mathrm{sw}}>0$。第一程以斜率 $v$ 扫至转折电位，第二程以斜率 $-v$ 扫回。FTacV 在两段电位上叠加同一个 $E_{\mathrm m}\sin(2\pi f_0t)$，其中 $E_{\mathrm m}$ 为峰值振幅，$f_0$ 为基频。两段沿同一个时钟计算正弦相位。[源码](https://github.com/ArsLinear/fdm-cvsim/blob/9e4f6b8755270f741ec8985a9a6233e1cde551ea/waveform.py)

| 函数 | 输入与输出 | 实现与调用 |
| --- | --- | --- |
| `cv_waveform_generator(initial, switch, scan_rate, dt)` | 起始和转折电位（V）、扫速（V/s）、正时间步长（s）→ 一维 CV 电位数组（V）。 | 用 `np.arange` 分别生成从 0 和从 $t_{\mathrm{sw}}$ 开始的两段时间，代入两段线性电位，再用 `np.concatenate` 拼接。 |
| `ftacv_waveform_generator(initial, switch, scan_rate, amplitude, frequency, dt)` | CV 参数、峰值振幅（V）、基频（Hz）→ 一维 FTacV 电位数组（V）。 | 分别计算两段直流扫描和对应的正弦扰动，逐点相加并拼接。 |

两函数只返回电位序列。若 $t_{\mathrm{sw}}/\Delta t$ 不是整数，两段拼接处相邻样点的间隔可能小于 $\Delta t$。将波形交给等时间步长的模拟器时，宜让转折点落在采样网格上。示例脚本直接从 CSV 取得电位。

## 正向模拟器

### `simulator/__init__.py`

包入口让调用者直接使用模型和模拟函数。网格接口从 `simulator.mesh` 导入。[源码](https://github.com/ArsLinear/fdm-cvsim/blob/9e4f6b8755270f741ec8985a9a6233e1cde551ea/simulator/__init__.py)

| 接口 | 输入与输出 | 实现与调用 |
| --- | --- | --- |
| `import simulator as sim` | 提供 `sim.OuterSphere`、`sim.faraday_sim`、`sim.edl_sim` 和 `sim.current_sim`。 | 从 `model.outer_sphere` 导出模型类，从 `simulator.py` 导出三个模拟函数。 |

### `simulator/mesh.py`

空间网格从 $x_0=0$ 开始。设 $L$ 为希望覆盖的溶液长度，$M$ 为区间数。等距网格令每个区间宽度为 $\Delta x$，取 $M=\lceil L/\Delta x\rceil$。指数网格令第 $i$ 格宽度为 $\Delta x_0r^i$，其中 $\Delta x_0$ 为第一格宽度，$r>1$ 为相邻格宽之比。程序取使格宽之和达到 $L$ 的最小整数 $M$，再累加各格宽得到节点。因此最后一个节点可能略超过 $L$。[源码](https://github.com/ArsLinear/fdm-cvsim/blob/9e4f6b8755270f741ec8985a9a6233e1cde551ea/simulator/mesh.py)

| 函数 | 输入与输出 | 实现与调用 |
| --- | --- | --- |
| `exponential_grid(length, dx0, ratio)` | 正目标长度 $L$（cm）、正首格宽 $\Delta x_0$（cm）、增长率 $r>1$ → 形状为 `(M+1,)` 的坐标数组（cm）。 | 用 `np.log1p` 根据等比数列求和式计算 $M=\left\lceil\log[1+(r-1)L/\Delta x_0]/\log r\right\rceil$，再累计 $M$ 个格宽并在首位补零。若需 $r=1$，调用 `uniform_grid`。 |
| `uniform_grid(length, dx)` | 正目标长度 $L$（cm）、正格宽 $\Delta x$（cm）→ 形状为 `(M+1,)` 的坐标数组（cm）。 | 将区间数向上取整，返回 $x_i=i\Delta x$，其中 $i=0,\ldots,M$。 |

### `simulator/model/outer_sphere.py`

考虑单步外球反应 $\mathrm{Ox}+ne^-\rightleftharpoons\mathrm{Red}$。$n$ 是电子转移数，$k_0$ 是标准速率常数，$\alpha$ 是电荷转移系数，$E^{0'}$ 是形式电位。设 $\mathbf c(x,t)=(c_{\mathrm{ox}}(x,t),c_{\mathrm{red}}(x,t))^{\mathsf T}$ 为浓度向量，$\mathbf c_{\mathrm s}=\mathbf c(0,t)$ 为电极表面浓度，$\mathbf D=\operatorname{diag}(D_{\mathrm{ox}},D_{\mathrm{red}})$ 为两种物种的扩散系数矩阵。定义过电位 $\eta=E-E^{0'}$ 和 $\beta=nF/(RT)$。源码固定气体常数 $R=8.314\ \mathrm{J\,mol^{-1}\,K^{-1}}$、法拉第常数 $F=96485\ \mathrm{C\,mol^{-1}}$、温度 $T=298.15\ \mathrm K$。

给定电位后，Butler–Volmer 净反应速率对两种表面浓度为线性。将反应系数写成矩阵 $\mathbf B(E)$，则表面边界和法拉第电流密度满足

$$
\begin{aligned}
\mathbf B(E)&=k_0
\begin{pmatrix}1\\-1\end{pmatrix}
\begin{pmatrix}\exp(-\alpha\beta\eta)&-\exp((1-\alpha)\beta\eta)\end{pmatrix},\\
\left.\mathbf D\,\partial_x\mathbf c\right|_{x=0}
&=\mathbf B(E)\mathbf c_{\mathrm s},\\
j_{\mathrm F}&=-nF[\mathbf B(E)\mathbf c_{\mathrm s}]_{\mathrm{ox}}.
\end{aligned}
$$

$[\mathbf B(E)\mathbf c_{\mathrm s}]_{\mathrm{ox}}>0$ 表示氧化态物种被净还原，因而 $j_{\mathrm F}<0$。模型对象的体相浓度同时用于整个浓度场的初值和模拟域远端的固定值。[源码](https://github.com/ArsLinear/fdm-cvsim/blob/9e4f6b8755270f741ec8985a9a6233e1cde551ea/simulator/model/outer_sphere.py)

| 函数或接口 | 输入与输出 | 实现与调用 |
| --- | --- | --- |
| `OuterSphere.__init__(self, waveform, k0, alpha, E0, n, D_ox, D_red, c_ox, c_red, C_dl)` | `waveform` 为界面电位序列，`E0` 对应 $E^{0'}$，`c_ox`、`c_red` 为体相浓度，`C_dl` 为面电容。调用 `OuterSphere(...)` 得模型对象。 | 构造函数将输入保存为同名属性。`faraday_sim` 读取扩散、界面反应和体相条件，`edl_current` 读取电位和面电容。 |
| `OuterSphere.diffusion_matrix(self)` | 模型中的 `D_ox`、`D_red` → 形状为 `(2, 2)` 的数组。 | 将两种扩散系数放在对角线上，供空间差分组装使用。 |
| `OuterSphere.interfacial_reaction_rate_coefficient_matrix(self, potential)` | 单个电位（V）→ 形状为 `(2, 2)` 的 $\mathbf B(E)$（cm/s）。 | 计算两个电位指数项，再用 `np.outer` 与化学计量向量 $(1,-1)^{\mathsf T}$ 作外积并乘以 $k_0$。 |
| `OuterSphere.butler_volmer(self, potential, concentration_surface)` | 单个电位（V）、按 `[ox, red]` 排列的表面浓度向量 → 标量 $j_{\mathrm F}$（A/cm²）。 | 调用速率矩阵方法，求 $\mathbf B(E)\mathbf c_{\mathrm s}$，取第一项乘以 $-nF$。 |

### `simulator/model/edl.py`

设 $C_{\mathrm{dl}}$ 为单位面积双电层电容。理想双电层的充放电电流密度为 $j_{\mathrm{dl}}=C_{\mathrm{dl}}\,\mathrm dE/\mathrm dt$。程序以等间隔电位差分逼近导数。内部样点用中心差分，首尾用单边差分，因此至少需要两个样点。[源码](https://github.com/ArsLinear/fdm-cvsim/blob/9e4f6b8755270f741ec8985a9a6233e1cde551ea/simulator/model/edl.py)

| 函数 | 输入与输出 | 实现与调用 |
| --- | --- | --- |
| `edl_current(model, dt)` | 含 `C_dl` 和 `waveform` 的模型、正时间步长（s）→ 形状为 `(N_t,)` 的双电层电流密度数组。 | 调用 `np.gradient(model.waveform, dt)` 求电位变化率，再乘以 `model.C_dl`。 |

### `simulator/assembly.py`

令 $\mathbf c_i^k$ 为时间点 $k$、空间节点 $i$ 处的浓度向量，$h_i=x_{i+1}-x_i$。对溶液内部的扩散方程 $\partial_t\mathbf c=\mathbf D\,\partial_x^2\mathbf c$，时间用后向欧拉差分，空间用非均匀三点差分。定义

$$
\mathbf A_i^-=\frac{2\Delta t\,\mathbf D}{h_{i-1}(h_{i-1}+h_i)},\qquad
\mathbf A_i^+=\frac{2\Delta t\,\mathbf D}{h_i(h_{i-1}+h_i)}.
$$

记 $\mathbf I$ 为 $2\times 2$ 单位矩阵。在内部节点 $1\le i\le N_x-2$，左、主、右三个系数块分别为 $-\mathbf A_i^-$、$\mathbf I+\mathbf A_i^-+\mathbf A_i^+$ 和 $-\mathbf A_i^+$，右端是上一时刻的 $\mathbf c_i^k$。在表面引入虚节点，用中心差分离散表面梯度，再代入反应边界消去虚节点。令 $\mathbf S=\Delta t\,\mathbf D/h_0^2$，$\mathbf Q_{k+1}=2\Delta t\,\mathbf B(E_{k+1})/h_0$，则表面行满足

$$
(\mathbf I+2\mathbf S+\mathbf Q_{k+1})\mathbf c_0^{k+1}
-2\mathbf S\,\mathbf c_1^{k+1}=\mathbf c_0^k .
$$

远端使用单位系数块，把 $\mathbf c_{N_x-1}^{k+1}$ 固定为体相浓度。未知量按空间节点排列，每个节点先放氧化态、后放还原态。整个方程组有 $2N_x$ 个未知量，随电位变化的是表面行的主对角块。[源码](https://github.com/ArsLinear/fdm-cvsim/blob/9e4f6b8755270f741ec8985a9a6233e1cde551ea/simulator/assembly.py)

| 函数 | 输入与输出 | 实现与调用 |
| --- | --- | --- |
| `boundary_matrix(interfacial_reaction_rate_coefficient_matrix, dx0, dt)` | 当前电位的 $\mathbf B(E)$、第一格宽度 $h_0$、步长 $\Delta t$ → 形状为 `(2, 2)` 的 $\mathbf Q$。 | 将反应系数矩阵乘以 $2\Delta t/h_0$，供 `big_matrix_refresh` 使用。`dx0` 取实际网格的第一格宽度。 |
| `concentration_matrix(initial_concentration, grid)` | 长度为 2 的体相浓度、含 $N_x$ 个节点的网格 → 形状为 `(N_x, 2)` 的初始浓度。 | 用 `np.tile` 把体相浓度复制到每个节点。 |
| `coefficient_matrix(diffusion_matrix, grid, dt)` | $\mathbf D$、严格递增的空间网格、正步长 → `(coefficient, surface_coefficient)`。 | 用 `np.diff(grid)` 求格宽。`coefficient` 形状为 `(2, N_x-2, 2, 2)`，首轴依次对应左、右系数。`surface_coefficient` 为形状 `(2, 2)` 的 $\mathbf S$。 |
| `big_matrix_init(coefficient, surface_coefficient)` | 内部系数块与表面扩散系数块 → 形状为 `(2*N_x, 2*N_x)` 的 LIL 稀疏矩阵。 | 填入内部三块、表面扩散块和远端单位块。表面反应块由 `big_matrix_refresh` 加入。 |
| `big_matrix_refresh(big_matrix, coefficient, flux)` | 现有稀疏矩阵、表面扩散块 $\mathbf S$、反应块 $\mathbf Q$ → 同一个稀疏矩阵。 | 就地覆盖左上角 $2\times2$ 块为 $\mathbf I+2\mathbf S+\mathbf Q$。这里的 `coefficient` 指表面系数，`flux` 已乘上时间与网格的缩放因子。 |
| `rhs_assemble(concentration, far_boundary)` | 上一步的 `(N_x, 2)` 浓度数组、长度为 2 的体相浓度 → 长度为 `2*N_x` 的右端向量。 | 复制旧浓度，把最后一行改为体相浓度，再按节点展平。 |

### `simulator/solver.py`

设 $\mathbf A$ 为当前时间步组装好的稀疏矩阵，$\mathbf b$ 为右端项。每次求解得到所有空间节点在新时刻的浓度向量 $\mathbf u$，满足 $\mathbf A\mathbf u=\mathbf b$。[源码](https://github.com/ArsLinear/fdm-cvsim/blob/9e4f6b8755270f741ec8985a9a6233e1cde551ea/simulator/solver.py)

| 函数 | 输入与输出 | 实现与调用 |
| --- | --- | --- |
| `solver(big_matrix, rhs)` | 稀疏系数矩阵、与矩阵行数一致的一维右端项 → 一维浓度解。 | 直接调用 `scipy.sparse.linalg.spsolve(big_matrix, rhs)`。调用者负责组装边界条件并还原数组形状。 |

### `simulator/simulator.py`

法拉第支路先初始化浓度场和稀疏矩阵。每个时间步用当前电位更新表面行，用上一步浓度组装右端，再解出当前浓度场，并用新表面浓度算 $j_{\mathrm F}$。双电层支路直接从电位序列得到 $j_{\mathrm{dl}}$。最后将两项相加得到总电流密度。[源码](https://github.com/ArsLinear/fdm-cvsim/blob/9e4f6b8755270f741ec8985a9a6233e1cde551ea/simulator/simulator.py)

| 函数 | 输入与输出 | 实现与调用 |
| --- | --- | --- |
| `faraday_sim(model, grid, dt)` | 模型、从 0 开始且至少含两个递增节点的网格、正时间步长 → `(time_array, faraday_current)`。两项形状均为 `(N_t,)`。 | 从模型取 $\mathbf D$、体相浓度和反应系数。初始化矩阵后转为 CSC 格式。从索引 1 开始调用 `boundary_matrix`、`big_matrix_refresh`、`rhs_assemble` 和 `solver`，将解还原成 `(N_x, 2)`，再调用 `model.butler_volmer`。`faraday_current[0]` 保留初始值 0。 |
| `edl_sim(model, dt)` | 模型、正时间步长 → 形状为 `(N_t,)` 的双电层电流密度。 | 将模型与步长交给 `edl_current`。 |
| `current_sim(model, grid, dt)` | 模型、网格、正时间步长 → `(time_array, total_current_density)`。模型电位至少有两个样点。 | 依次调用 `faraday_sim` 和 `edl_sim`，返回 $t_k=k\Delta t$ 和 $j_{\mathrm F}+j_{\mathrm{dl}}$。首个总电流密度样点来自双电层项。 |

## 导出与网格灵敏度

### `plot_ftacv.py` 与 `plot_voltammetry.py`

这两个脚本将同一条模拟记录整理成两种表格。一种保留时间、电位、电流密度三列，便于继续做信号处理。另一种保留电位与电流密度两列，便于画伏安曲线。两个脚本导入本地 `ftacvsim.py` 以获得数组，执行对应的 `main()` 后把 CSV 写到当前工作目录。[plot_ftacv.py](https://github.com/ArsLinear/fdm-cvsim/blob/9e4f6b8755270f741ec8985a9a6233e1cde551ea/plot_ftacv.py) · [plot_voltammetry.py](https://github.com/ArsLinear/fdm-cvsim/blob/9e4f6b8755270f741ec8985a9a6233e1cde551ea/plot_voltammetry.py)

| 函数 | 输入与输出 | 实现与调用 |
| --- | --- | --- |
| `plot_ftacv.main()` | 无参数。读取同长的 `time_array`、`waveform_array` 和 `current`，写出 `ftacvsim.csv`，返回 `None`。 | 用 `zip(..., strict=True)` 配对样点，用 `csv.writer` 写出时间（s）、电位（V）和电流密度（A/cm²），然后打印行数。 |
| `plot_voltammetry.main()` | 无参数。读取同长的 `waveform_array` 与 `current`，写出 `voltammogram.csv`，返回 `None`。 | 同样按索引配对样点，写出电位（V）和电流密度（A/cm²），然后打印行数。 |

两个脚本都用 `Path("./result").with_name(...)` 生成路径。`with_name` 把末段 `result` 替换为文件名，因此目标文件位于当前目录。直接执行脚本会调用 `main()`，并覆盖同名 CSV。

### `sensitivity.py`

本文件估计第一格宽度引起的离散误差有多大可能表现为参数偏移。设 $\mathbf j_{\mathrm{ref}}$ 为细网格上的总电流密度序列，$\mathbf j_{\mathrm{grid}}$ 为候选网格上的序列，基准参数为 $\boldsymbol\theta=(k_0,\alpha,E^{0'},C_{\mathrm{dl}})$。设 $\delta_j$ 为第 $j$ 个参数的相对变化，即扰动后参数为 $\theta_j(1+\delta_j)$。程序在参考网格上每次只改变一个参数，用 $\theta_j(1\pm\epsilon)$ 各模拟一次，其中 $\epsilon=10^{-4}$。两条电流密度序列之差除以 $2\epsilon$，得到雅可比矩阵 $\mathbf J$ 的第 $j$ 列。若基准参数为零，相对扰动对应的列也为零。

在局部线性近似下，用参数扰动产生的电流变化 $\mathbf J\boldsymbol\delta$ 抵消候选网格误差。程序用普通最小二乘求

$$
\widehat{\boldsymbol\delta}=\arg\min_{\boldsymbol\delta}
\left\|\mathbf J\boldsymbol\delta+(\mathbf j_{\mathrm{grid}}-\mathbf j_{\mathrm{ref}})\right\|_2^2.
$$

这一步使用模拟曲线比较网格离散误差，所得是相对参数偏移的线性估计。脚本取 $10^{-10}$ cm 为参考首格宽度，遍历 $10^{-3-0.25k}$ cm（$k=0,\ldots,15$）作为候选宽度。四项相对偏移的绝对值均不超过 0.5% 时打印 `PASS`。参考电流、雅可比和扫描循环均在模块顶层执行。[源码](https://github.com/ArsLinear/fdm-cvsim/blob/9e4f6b8755270f741ec8985a9a6233e1cde551ea/sensitivity.py)

| 函数 | 输入与输出 | 实现与调用 |
| --- | --- | --- |
| `grid(dx0)` | 正首格宽度（cm）→ 一维坐标数组（cm）。 | 以 `base.length` 为目标长度，调用 `exponential_grid`。保留严格小于目标长度的节点，再追加精确的远端坐标，使不同候选网格共享同一个远端边界。 |
| `sim_current(dx0, model)` | 首格宽度、模型 → 形状为 `(N_t,)` 的总电流密度数组。 | 调用 `grid(dx0)`，再用从基准脚本取得的 `DT` 调用 `sim.current_sim`，取电流密度数组。 |
| `model_with_parameter(name, value)` | 参数属性名、新值（沿用原参数单位）→ 修改该属性后的模型副本。 | 对 `base.model` 做浅复制，用 `setattr` 设置一个属性。电位数组等其他对象仍由副本共享。 |
| `calculate_jacobian(dx0_ref)` | 参考首格宽度（cm）→ `(reference_current, jacobian)`，形状分别为 `(N_t,)` 和 `(N_t, 4)`。 | 先模拟基准曲线，再按 `k0`、`alpha`、`E0`、`C_dl` 的顺序逐项做相对扰动。每次由 `model_with_parameter` 生成模型，再用 `sim_current` 计算曲线并作中心差分。 |
| `estimate_parameter_shift(dx0, reference_current, jacobian)` | 候选首格宽度、参考电流密度数组、参考雅可比矩阵 → 形状为 `(4,)` 的无量纲相对偏移。 | 用基准模型调用 `sim_current` 求候选曲线，以参考曲线减候选曲线为右端，调用 `np.linalg.lstsq` 求解。 |
