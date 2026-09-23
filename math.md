# 数学推理

## 理想外球反应的模型建立

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

1. 电极反应的消耗或生成。这一部分仅发生在电极表面，可以通过法拉第定律来描述，即：
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
因此，求循环伏安的伏安曲线，本质上就是解这个偏微分方程组。

## 有限差分法（FDM）求解方程组

### 有限差分法（FDM）简介

给定定义域在$[0,L]$上的函数$f(x)$，将定义域分割为若干个长度为$\Delta x$的网格。第$i$个网格前后端点分别为$f_{i-1}(x)$与$f_i(x)$，则：
$$
\dv{f}{x} \approx \frac{f_{i+1}(x)-f_{i-1}(x)}{2\Delta x}\\\\
\dv[2]{f}{x} \approx \frac{f_{i+1}(x)-2f_{i}(x)+f_{i-1}(x)}{\Delta x^2}
$$
所以，将Fick第二定律在时间和空间上离散就可以得到：
$$
\frac{\mathbf{c}_{i}^{j+1}-\mathbf{c}_{i}^{j}}{\Delta t}

=

\mathbf{D}\frac{\mathbf{c}_{i+1}^j-2\mathbf{c}_{i}^j+\mathbf{c}_{i-1}^j}{\Delta x^2}
$$
因此，我们可以通过这一时刻的浓度空间分布数值求解出下一时刻的浓度空间分布：
$$
\mathbf{c}_{i}^{j+1}=\frac{\mathbf{D}\Delta t}{\Delta x^2}(\mathbf{c}_{i+1}^j-2\mathbf{c}_{i}^j+\mathbf{c}_{i-1}^j)+\mathbf{c}_{i}^{j}
$$
这就是FDM的显式推进。

显式推进认为每一个时步后的状态是由扩散方程从前一个时步推断出来的，而另一种思路则认为，扩散方程在每一个时间点中对全体均成立，也就是说，后一个时刻的浓度分布相互耦合并且应该由前一个时刻浓度分布全体决定。令$\mathbf{R}=\frac{\mathbf{D}\Delta t}{\Delta x^2}$：
$$
\mathbf{c}_{i}^{j+1}-\mathbf{c}_{i}^{j}=\mathbf{R}(\mathbf{c}_{i+1}^{j+1}-2\mathbf{c}_{i}^{j+1}+\mathbf{c}_{i-1}^{j+1})
$$
所以
$$
-\mathbf{R}\mathbf{c}_{i-1}^{j+1} + (\mathbf{I} + 2\mathbf{R})\mathbf{c}_{i}^{j+1} - \mathbf{R}\mathbf{c}_{i+1}^{j+1} = \mathbf{c}_{i}^{j}
$$
我们联立所有的$i$，解方程组：
$$
\begin{pmatrix}
\mathbf{I} + 2\mathbf{R} & -2\mathbf{R} & 0 & \cdots \\
-\mathbf{R} & \mathbf{I} + 2\mathbf{R} & -\mathbf{R} & \cdots \\
0 & -\mathbf{R} & \mathbf{I} + 2\mathbf{R} & \ddots \\
\vdots & \vdots & \ddots & \ddots
\end{pmatrix}
\begin{pmatrix}
\mathbf{c}_0^{j+1} \\
\mathbf{c}_1^{j+1} \\
\mathbf{c}_2^{j+1} \\
\vdots
\end{pmatrix}
=
\begin{pmatrix}
\mathbf{c}_0^j \\
\mathbf{c}_1^j \\
\mathbf{c}_2^j \\
\vdots
\end{pmatrix}
$$
则可以解得$j+1$时的物质浓度分布。这就是FDM的隐式推进。

### 隐式推进解方程组

上述计算方法是针对无源项、无流项的Fick扩散过程。如果界面有浓度流（比如说我们的法拉第过程），该如何使用FDM求解呢？在这个过程中，电极表面$x=0$处有：
$$
\mathbf{D}\pdv{\mathbf{c}}{x}
=k_0\mathbf{s}\mathbf{A}^{\mathsf T}\mathbf{c}
$$
对空间离散可以得到：
$$
\mathbf{D}\frac{\mathbf{c}_{1}^{j+1}-\mathbf{c}_{-1}^{j+1}}{2\Delta x}
=k_0\mathbf{s}\mathbf{A}_{j+1}^{\mathsf T}\mathbf{c}_{0}^{j+1}
$$
在电极界面上的Fick扩散定律，有
$$
\mathbf{c}_{0}^{j+1}-\mathbf{c}_{0}^{j}=\mathbf{R}(\mathbf{c}_{1}^{j+1}-2\mathbf{c}_{0}^{j+1}+\mathbf{c}_{-1}^{j+1})
$$
联立(16)(17)得：
$$
\mathbf{c}_{0}^{j+1}-\mathbf{c}_{0}^{j}=2\mathbf{R}(\mathbf{c}_{1}^{j+1}-2\mathbf{c}_{0}^{j+1})+\frac{2k_0\Delta t}{\Delta x} \mathbf{s}\mathbf{A}_{j+1}^{\mathrm{T}}\mathbf{c}_{-1}^{j+1}
$$
令$\mathbf{Q}_{j+1} = \frac{2k_0\Delta t}{\Delta x} \mathbf{s}\mathbf{A}_{j+1}^{\mathrm{T}}$则：
$$
(\mathbf{I} + 2\mathbf{R} - \mathbf{Q}_{j+1})\mathbf{c}_0^{j+1} - 2\mathbf{R}\mathbf{c}_1^{j+1} = \mathbf{c}_0^j
$$
所以我们需要求解的方程组变成了：
$$
\begin{pmatrix}
\mathbf{I} + 2\mathbf{R} - \mathbf{Q} & -2\mathbf{R} & 0 & \cdots \\
-\mathbf{R} & \mathbf{I} + 2\mathbf{R} & -\mathbf{R} & \cdots \\
0 & -\mathbf{R} & \mathbf{I} + 2\mathbf{R} & \ddots \\
\vdots & \vdots & \ddots & \ddots
\end{pmatrix}
\begin{pmatrix}
\mathbf{c}_0^{j+1} \\
\mathbf{c}_1^{j+1} \\
\mathbf{c}_2^{j+1} \\
\vdots
\end{pmatrix}
=
\begin{pmatrix}
\mathbf{c}_0^j \\
\mathbf{c}_1^j \\
\mathbf{c}_2^j \\
\vdots
\end{pmatrix}
$$
