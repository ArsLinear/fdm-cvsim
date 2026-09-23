# 关于 `math.md` 表面边界条件矩阵符号的修正

## 结论

表面行的正确形式是

$$
(\mathbf{I}+2\mathbf{R}+\mathbf{Q})\,\mathbf{c}_0^{j+1}-2\mathbf{R}\,\mathbf{c}_1^{j+1}=\mathbf{c}_0^{j},
\qquad
\mathbf{Q}=\frac{2k_0\Delta t}{\Delta x}\,\mathbf{s}\mathbf{A}^{\mathsf T},
$$

即应是 $+\mathbf{Q}$，且 $\mathbf{s}\mathbf{A}^{\mathsf T}$ 的对角元为正。`math.md` 最终矩阵给出的 $\mathbf{I}+2\mathbf{R}-\mathbf{Q}$ 符号相反。

## 推导

由 `math.md` (16)(17)：

$$
\mathbf{D}\frac{\mathbf{c}_1^{j+1}-\mathbf{c}_{-1}^{j+1}}{2\Delta x}=k_0\mathbf{s}\mathbf{A}_{j+1}^{\mathsf T}\mathbf{c}_0^{j+1}
\tag{16}
$$

$$
\mathbf{c}_0^{j+1}-\mathbf{c}_0^{j}
=\mathbf{R}\,(\mathbf{c}_1^{j+1}-2\mathbf{c}_0^{j+1}+\mathbf{c}_{-1}^{j+1})
\tag{17}
$$

式 (16) 解出幽灵点：

$$
\mathbf{D}\,\mathbf{c}_{-1}^{j+1}
=\mathbf{D}\,\mathbf{c}_1^{j+1}-2\Delta x\,k_0\,\mathbf{s}\mathbf{A}^{\mathsf T}\mathbf{c}_0^{j+1}
\quad\Longrightarrow\quad
\mathbf{c}_{-1}^{j+1}
=\mathbf{c}_1^{j+1}-2\Delta x\,\mathbf{D}^{-1}k_0\,\mathbf{s}\mathbf{A}^{\mathsf T}\mathbf{c}_0^{j+1}
$$

代入 (17)：

$$
\mathbf{c}_0^{j+1}-\mathbf{c}_0^{j}
=\mathbf{R}\left(2\mathbf{c}_1^{j+1}-2\mathbf{c}_0^{j+1}-2\Delta x\,\mathbf{D}^{-1}k_0\,\mathbf{s}\mathbf{A}^{\mathsf T}\mathbf{c}_0^{j+1}\right)
$$

整理，并注意 $\mathbf{R}\mathbf{D}^{-1}\Delta x=\dfrac{\Delta t}{\Delta x}\mathbf{I}$：

$$
\bigl(\mathbf{I}+2\mathbf{R}+\underbrace{\tfrac{2k_0\Delta t}{\Delta x}\mathbf{s}\mathbf{A}^{\mathsf T}}_{\mathbf{Q}}\bigr)\mathbf{c}_0^{j+1}-2\mathbf{R}\,\mathbf{c}_1^{j+1}=\mathbf{c}_0^{j}
$$

## `math.md` 的两处笔误

1. 式 (18) 右端最后一项写成了 $\dfrac{2k_0\Delta t}{\Delta x}\mathbf{s}\mathbf{A}^{\mathsf T}\mathbf{c}_{-1}^{j+1}$，应当是 $\mathbf{c}_0^{j+1}$（幽灵点已在上一步消去）。
2. 最终矩阵把这一项写成 $-\mathbf{Q}$，应为 $+\mathbf{Q}$。

## 为什么必须是 $+\mathbf{Q}$

- **M-矩阵 / 非负性**：$\mathbf{Q}$ 的两个对角元为 $+$（$\left[\mathbf{s}\mathbf{A}^{\mathsf T}\right]_{11}=a_\text{ox}>0$，$\left[\mathbf{s}\mathbf{A}^{\mathsf T}\right]_{22}=a_\text{red}>0$），使 $\mathbf{I}+2\mathbf{R}+\mathbf{Q}$ 保持对角占优；若写成 $-\mathbf{Q}$，对角被削弱甚至反号，解出的浓度会出现负值。
- **Nernst 极限自洽**：$Q\to\infty$ 时（无限快动力学）表面行由 $\mathbf{Q}\mathbf{c}_0\approx\mathbf{0}$ 主导，其零空间给出 $c_\text{ox}/c_\text{red}=e^{\xi}$（$\xi=\frac{nF}{RT}(E-E^{0'})$），正是 Nernst 方程；若符号反了，该极限会变成浓度发散的双曲解而非平衡态。

## 对代码的影响

`cvsim.py` 中 `M[0:2, 0:2] = Mii + (2*k0*dt/dx) * np.outer([1,-1], [a_ox, -a_red])` 即为 $(\mathbf{I}+2\mathbf{R}+\mathbf{Q})$ 的左上块。
