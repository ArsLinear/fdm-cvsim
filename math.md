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

c_\text{ox}^\text{s}\exp\left(-\frac{\alpha nF}{RT}\right)(E-E^{0'})

-

c_\text{red}^\text{s}\exp\left(\frac{(1-\alpha) nF}{RT}\right)(E-E^{0'})

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
D_\text{ox} \\
D_\text{red}
\end{pmatrix}

\\\\

\mathbf{A} 
=
\begin{pmatrix}
\exp[\left(-\frac{\alpha nF}{RT}\right)(E-E^{0'})] \\
\exp[\left(\frac{(1-\alpha) nF}{RT}\right)(E-E^{0'})]
\end{pmatrix}\\
$$
那么上述边界条件可以写成线性形式：
$$
\left\{

\begin{aligned}

\frac{\partial \mathbf{c}}{\partial t}
&=-\frac{\mathbf{A}\cdot\mathbf{c}}{nF},
&& x=0 \\[6pt]

\frac{\partial \mathbf{c}}{\partial t}
&=\mathbf{D}^{\mathsf T}\nabla^2\mathbf{c},
&& 0<x<L \\[6pt]

\mathbf{c}
&=\mathbf{c}_{\mathrm{bulk}},
&& x=L

\end{aligned}

\right.
$$
因此，求循环伏安的伏安曲线，本质上就是解这个偏微分方程组。
