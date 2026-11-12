# 数学实现

## 理想双电层电容的充放电电流

理想外球反应的非法拉第过程可以归属为电极表面双电层的充放电。理想双电层可以视为一个理想电容器，与Butler-Volmer元件并联，因此非法拉第电流可以计算为：
$$
i=C_\text{dl}\dv{E}{t}
$$

## 理想外球反应的法拉第电流

循环伏安过程中，一段直流扫描的激励信号可以写为：
$$
E(t)=E_0+vt
$$
对于外球反应，我们先考虑电荷转移过程。若极化均用于电荷转移，则根据Butler-Vomle方程：
$$
j=nk_0F

\left[

c_\text{ox}^\text{s}\exp\left(-\frac{\alpha nF}{RT}(E-E^{0'})\right)

-

c_\text{red}^\text{s}\exp\left(\frac{(1-\alpha) nF}{RT}(E-E^{0'})\right)

\right]
$$
将激励信号公式代入即可。

然而在这个体系下，$c^\text{s}$在随着电势的变化而不断变化。所以若想求循环伏安的解析表达，我们需要计算表面物种浓度随电势的变化。

表面物种浓度受到两个因素的影响：

1. 电极反应的消耗或生成。法拉第过程生成或消耗的浓度流等于界面浓度流。这一部分仅发生在电极表面，可以通过法拉第定律来描述，即：
   $$
   \pdv{c^\text{s}}{t}=-\frac{j}{nF}
   $$

2. 与体相溶液之间的扩散。这一部分发生在电极附近一段距离，可以通过Fick第二定律来描述，即：
   $$
   \pdv{c}{t}=D\nabla^2c
   $$

如果我们令：
$$
\mathbf{c} 
=
\begin{pmatrix}
c_\text{ox} \\
c_\text{red}
\end{pmatrix}

\\\\

\mathbf{D} 
=
\begin{pmatrix}
D_\text{ox}&0 \\
0&D_\text{red}
\end{pmatrix}

\\\\

\mathbf{A} 
=
\begin{pmatrix}
\exp[\left(-\frac{\alpha nF}{RT}\right)(E-E^{0'})] \\
-\exp[\left(\frac{(1-\alpha) nF}{RT}\right)(E-E^{0'})]
\end{pmatrix}

\\\\

\mathbf{s}
=
\begin{pmatrix}
1\\
-1
\end{pmatrix}\\
$$
那么上述边界条件可以写成线性形式：
$$
\left\{

\begin{aligned}

\mathbf{D}\pdv{\mathbf{c}}{x}
&=k_0\mathbf{s}\mathbf{A}^{\mathsf T}\mathbf{c},
&& x=0 \\[6pt]

\frac{\partial \mathbf{c}}{\partial t}
&=\mathbf{D}\pdv[2]{\mathbf{c}}{x},
&& 0<x<L \\[6pt]

\mathbf{c}
&=\mathbf{c}_{\mathrm{bulk}},
&& x=L

\end{aligned}

\right.
$$
我们将电极表面的网格节点定义为$x_0$，电极表面节点到体相第一个节点的距离定义为$\Delta x_0$，网格节点数量为$N$，那么最后一个网格距离就是$\Delta x_{N-1}$。指数网格我们采用：
$$
\Delta x_{i+1} =\Delta x_0\cdot q^n
$$
所以对于网格的划分，我们只需要确定两个参数：$\Delta x_0$与$q$。这两个量具有明确的物理含义：

- 电极表面节点到体相第一个节点的距离$\Delta x_0$，其物理含义是电极附近的分辨率；
- 相邻间距的比例$q$，其物理含义是离开电极后网格变疏的速度。

指数网格下二阶导数的离散会发生一些变化。将函数在$x_{i}$附近泰勒展开可以得到：
$$
\begin{aligned}
f_{i+1} &= f_i + \Delta x_i f_i' + \frac{\Delta x_i^2}{2} f_i'' + \frac{\Delta x_i^3}{6} f_i''' + \cdots  \\

f_{i-1} &= f_i - \Delta x_{i-1} f_i' + \frac{\Delta x_{i-1}^2}{2} f_i'' - \frac{\Delta x_{i-1}^3}{6} f_i''' + \cdots 
\end{aligned}
$$
忽略三阶以上无穷小量，我们能通过上面两个式子将一阶导数消去得到：
$$
\dv[2]{f_{i+1}}{x}
\approx
\frac{2}{\Delta x_{i-1} + \Delta x_i}
\left[
  \frac{f_{i+1} - f_i}{\Delta x_i}
  -
  \frac{f_i - f_{i-1}}{\Delta x_{i-1}}
\right]
$$
那么这个就是线性网格下的离散处理方法。这个式子有着明确的物理意义。$({f_{i+1} - f_i})/{\Delta x_i}$是一阶导数，而$(\Delta x_{i-1} + \Delta x_i)/2$就是二者网格的中心距离。所以当$\Delta x_{i-1} =\Delta x_i$时，该离散形式就会直接退化为均匀网格的二阶导数离散形式。

对于体相的扩散控制方程，我们依然对时间均匀离散，而改变空间的离散方式，可以得到：
$$
\frac{\mathbf{c}_i^{j+1} - \mathbf{c}_i^j}{\Delta t}
=
\frac{2\mathbf{D}}{\Delta x_{i-1} + \Delta x_i}
\left[
  \frac{\mathbf{c}_{i+1}^{j+1} - \mathbf{c}_i^{j+1}}{\Delta x_i}
  -
  \frac{\mathbf{c}_i^{j+1} - \mathbf{c}_{i-1}^{j+1}}{\Delta x_{i-1}}
\right]
$$
我们定义常数：
$$
\begin{aligned}

\mathbf{R}_i^- = \frac{2 \mathbf{D} \Delta t}{\Delta x_{i-1} (\Delta x_{i-1} + \Delta x_i)} \\[6pt]

\mathbf{R}_i^+ = \frac{2 \mathbf{D} \Delta t}{\Delta x_i (\Delta x_{i-1} + \Delta x_i)}\

\end{aligned}
$$
这样所有常数就可以被统一归到系数矩阵中：
$$
\mathbf{c}_i^{j+1} - \mathbf{c}_i^j

=

\mathbf{R}_i^-(\mathbf{c}_{i+1}^{j+1} - \mathbf{c}_i^{j+1})

-

\mathbf{R}_i^+(\mathbf{c}_i^{j+1} - \mathbf{c}_{i-1}^{j+1})
$$
化简得到
$$
-\mathbf{R}_i^- \mathbf{c}_{i-1}^{j+1} + (\mathbf{I} + \mathbf{R}_i^- + \mathbf{R}_i^+) \mathbf{c}_i^{j+1} - \mathbf{R}_i^+ \mathbf{c}_{i+1}^{j+1} = \mathbf{c}_i^j
$$
现在来处理边界流。实际上，这个控制方程没有涉及到浓度函数的二阶导数，所以这个控制方程可以完全套用原来的结果。同时，我们将Ghost Point放在与第一个节点关于镜像对称的地方，也就是说，在第一个节点我们可以套用原有的体相扩散表达式。所以离散边界流的控制方程可以写作：
$$
\mathbf{D} \frac{\mathbf{c}_1^{j+1} - \mathbf{c}_{-1}^{j+1}}{2 \Delta x_0} = k_0 s \mathbf{A}_{j+1}^{\mathrm{T}} \mathbf{c}_0^{j+1}
$$
这样，电极表面的控制方程就可以完全的复用原有的结果，即：
$$
(\mathbf{I} + 2\mathbf{R}_0 + \mathbf{Q}_{j+1})\mathbf{c}_0^{j+1} - 2\mathbf{R}_0\mathbf{c}_1^{j+1} = \mathbf{c}_0^j
$$
其中，$\mathbf{Q}_{j+1} = \frac{2k_0\Delta t}{\Delta x_0} \mathbf{s}\mathbf{A}_{j+1}^{\mathrm{T}}$，$\mathbf{R}_0=\frac{\mathbf{D}\Delta t}{\Delta x_0^2}$。这样，整个矩阵方程就可以改写为：
$$
\begin{pmatrix}
\mathbf I+2\mathbf R_0+\mathbf Q_{j+1}
& -2\mathbf R_0 & 0 & \cdots & 0 \\
-\mathbf R_1^-
& \mathbf I+\mathbf R_1^-+\mathbf R_1^+
& -\mathbf R_1^+ & \cdots & 0 \\
0 & \ddots & \ddots & \ddots & \vdots \\
\vdots & \cdots & -\mathbf R_{N-1}^-
& \mathbf I+\mathbf R_{N-1}^-+\mathbf R_{N-1}^+
& -\mathbf R_{N-1}^+ \\
0 & \cdots & 0 & 0 & \mathbf I
\end{pmatrix}
\begin{pmatrix}
\mathbf c_0^{j+1} \\
\mathbf c_1^{j+1} \\
\vdots \\
\mathbf c_{N-1}^{j+1} \\
\mathbf c_N^{j+1}
\end{pmatrix}
=
\begin{pmatrix}
\mathbf c_0^j \\
\mathbf c_1^j \\
\vdots \\
\mathbf c_{N-1}^j \\
\mathbf c_{\mathrm{bulk}}
\end{pmatrix}
$$
