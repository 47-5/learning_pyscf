# Basis 中的数字如何进入 GTO 公式

这份笔记只讨论 orbital basis，也就是通常写在 `mol.basis` 里的 AO 基组。ECP 和 pseudo 的数字格式是另一套东西，不放在这里混讲。

## 1. 从分子轨道展开开始

量子化学程序通常把分子轨道写成 AO basis 的线性组合：

$$
\psi_i(\mathbf r) = \sum_\mu C_{\mu i}\,\chi_\mu(\mathbf r)
$$

这里：

- $\psi_i$ 是第 $i$ 个分子轨道。
- $\chi_\mu$ 是第 $\mu$ 个 AO basis function。
- $C_{\mu i}$ 是 MO coefficient，也就是 SCF 过程中求出来的未知量。

基组文件里的数字定义的是 $\chi_\mu$。也就是说，basis 不是直接给出分子轨道，而是给出一套可用于展开分子轨道的原子中心函数。

## 2. 一个 primitive Gaussian 的公式

PySCF 的常规 Gaussian basis 最终对应 Gaussian type orbital, 简写 GTO。以原子 $A$ 为中心，令

$$
\mathbf r_A = \mathbf r - \mathbf R_A, \qquad r_A = |\mathbf r_A|
$$

球谐形式下，一个 primitive GTO 可以理解为

$$
g_{p l m}^{A}(\mathbf r)
= N_l(\alpha_p)\, r_A^l e^{-\alpha_p r_A^2}\,Y_{lm}(\hat{\mathbf r}_A)
$$

其中：

- $p$ 是 primitive 编号。
- $l$ 是角动量量子数，$s,p,d,f,\cdots$ 分别对应 $l=0,1,2,3,\cdots$。
- $m$ 是同一个 $l$ 下的角向分量。球谐基下共有 $2l+1$ 个分量。
- $\alpha_p$ 是 Gaussian exponent，控制函数径向衰减快慢。
- $N_l(\alpha_p)$ 是 primitive normalization factor。

在 PySCF 里，primitive radial part 的归一化因子由 `gto_norm(l, alpha)` 计算：

$$
N_l(\alpha)
= \frac{1}{\sqrt{\int_0^\infty r^{2l+2}e^{-2\alpha r^2}\,dr}}
$$

源码里也等价写成：

$$
N_l(\alpha)
= \sqrt{\frac{2^{2l+3}(l+1)!(2\alpha)^{l+3/2}}{(2l+2)!\sqrt{\pi}}}
$$

要注意这里的 $N_l(\alpha)$ 是 radial part 的归一化。真实 AO 的完整归一化还和角向函数约定有关，PySCF/libcint 会按照 spherical 或 Cartesian 的积分约定处理。

## 3. Cartesian 形式怎么看

如果使用 Cartesian Gaussian，同一个 shell 的角向部分不是 $Y_{lm}$，而是 Cartesian polynomial：

$$
g_{pabc}^{A}(\mathbf r)
= N_{abc}(\alpha_p)\,x_A^a y_A^b z_A^c e^{-\alpha_p r_A^2}
$$

其中

$$
a+b+c=l
$$

因此同一个 $l$ shell 里 Cartesian 函数个数是

$$
\frac{(l+1)(l+2)}{2}
$$

而 spherical 函数个数是

$$
2l+1
$$

例如：

| shell | $l$ | spherical 个数 | Cartesian 个数 |
| --- | ---: | ---: | ---: |
| s | 0 | 1 | 1 |
| p | 1 | 3 | 3 |
| d | 2 | 5 | 6 |
| f | 3 | 7 | 10 |

这就是为什么同样写一个 d shell，有些程序会讨论 5D/6D，有些会讨论 spherical/cartesian。

## 4. contraction: 基组文件中的 coefficient 是什么

实际计算很少直接把每个 primitive 都当成独立 AO。通常会把多个 primitive 线性组合成 contracted Gaussian：

$$
\chi_{k l m}^{A}(\mathbf r)
= \mathcal N_{kl}^{\mathrm{contr}}
\sum_{p=1}^{n_{\mathrm{prim}}}
 c_{pk}\,N_l(\alpha_p)\,r_A^l e^{-\alpha_p r_A^2}Y_{lm}(\hat{\mathbf r}_A)
$$

这里：

- $n_{\mathrm{prim}}$ 是 primitive 数量。
- $k$ 是 contraction 编号。
- $c_{pk}$ 是基组表中看到的 contraction coefficient。
- $\mathcal N_{kl}^{\mathrm{contr}}$ 是 contracted AO 的整体归一化因子。PySCF 默认会做这一步，因为 `NORMALIZE_GTO = True`。

如果一个 shell 只有一列 coefficient，那么 $k=1$。这句话更准确地说是：这个 shell 只有一个 contracted radial function，而不是说它只产生一个 AO。可以先把径向部分单独抽出来：

$$
R_{k l}^{A}(r_A)
= \mathcal N_{kl}^{\mathrm{contr}}
\sum_{p=1}^{n_{\mathrm{prim}}} c_{pk}\,N_l(\alpha_p)\,r_A^l e^{-\alpha_p r_A^2}
$$

完整 AO 还要再乘上角向部分：

$$
\chi_{k l m}^{A}(\mathbf r)=R_{kl}^{A}(r_A)Y_{lm}(\hat{\mathbf r}_A)
$$

所以“几列 coefficient”和“几个角向分量”是两件事：

| 层次 | 由谁决定 | 数量 |
| --- | --- | --- |
| contracted radial function | coefficient 有几列，即 `nctr` | $n_{\mathrm{ctr}}$ |
| angular functions | shell 的 $l$ 和 spherical/cartesian 约定 | spherical: $2l+1$；Cartesian: $(l+1)(l+2)/2$ |
| 最终 AO 数 | 两者相乘 | `nctr * angular_count` |

例如一个 p shell，`l = 1` 且只有一列 coefficient。为了避免把“primitive 数量”和 $p_x,p_y,p_z$ 的三个方向混在一起，这里假设它有 $M$ 个 primitive，而不是恰好写成 3 个：

```python
[1,
  [alpha_1, c_1],
  [alpha_2, c_2],
  ...,
  [alpha_M, c_M]]
```

这里的 $M=n_{\mathrm{prim}}$，表示 primitive Gaussian 的个数。它和 p shell 的三个角向分量不是同一个概念。

它只有一个径向收缩函数：

$$
R_{1p}^{A}(r_A)
= \mathcal N^{\mathrm{contr}}
\sum_{q=1}^{M}c_qN_1(\alpha_q)r_Ae^{-\alpha_qr_A^2}
$$

这里用 $q$ 做 primitive 编号，是为了避免和 p shell 的字母 p 混淆。

但是它会产生三个 p 型 AO。用 Cartesian 语言写就是：

$$
\begin{aligned}
\chi_{p_x}^{A}(\mathbf r) &= R_{1p}^{A}(r_A)\,\frac{x_A}{r_A} \\
\chi_{p_y}^{A}(\mathbf r) &= R_{1p}^{A}(r_A)\,\frac{y_A}{r_A} \\
\chi_{p_z}^{A}(\mathbf r) &= R_{1p}^{A}(r_A)\,\frac{z_A}{r_A}
\end{aligned}
$$

也常等价地写成 Cartesian primitive polynomial 的形式：

$$
\chi_{p_x}^{A}(\mathbf r)
= \mathcal N^{\mathrm{contr}}
\sum_{q=1}^{M}c_qN_x(\alpha_q)x_Ae^{-\alpha_qr_A^2}
$$

$\chi_{p_y}$ 和 $\chi_{p_z}$ 只是把 $x_A$ 换成 $y_A,z_A$。如果使用 spherical AO，则不是直接叫 $p_x,p_y,p_z$，而是三个 real spherical harmonic 组合；但数量仍然是 3。

这组三个函数也可以更紧凑地写成向量形式：

$$
\boldsymbol\chi_{p}^{A}(\mathbf r)
=
\begin{bmatrix}
\chi_{p_x}^{A}(\mathbf r) \\
\chi_{p_y}^{A}(\mathbf r) \\
\chi_{p_z}^{A}(\mathbf r)
\end{bmatrix}
=
\mathcal N^{\mathrm{contr}}
\sum_{q=1}^{M}c_qN_1(\alpha_q)e^{-\alpha_qr_A^2}
\begin{bmatrix}
x_A \\
y_A \\
z_A
\end{bmatrix}
$$

这里最重要的观察是：同一个径向收缩系数 $c_q$ 同时作用在三个角向分量上。换句话说，p shell 的 contraction 不是分别给 $p_x,p_y,p_z$ 三套不同的 coefficient，而是同一组 coefficient 乘上一个三维角向向量。

为了和后面的矩阵 $\mathbf B$ 无缝衔接，可以把每个 primitive p shell 也写成一个三维列向量：

$$
\mathbf b_q^{A,p}(\mathbf r)
=
N_1(\alpha_q)e^{-\alpha_qr_A^2}
\begin{bmatrix}
x_A \\
y_A \\
z_A
\end{bmatrix},
\qquad q=1,2,\cdots,M
$$

所以 $\mathbf b_1,\mathbf b_2,\cdots,\mathbf b_M$ 的下标表示第几个 primitive Gaussian，而不是 $p_x,p_y,p_z$。每一个 $\mathbf b_q$ 本身都是一个三维列向量，已经包含三个方向。

于是单个 contraction 就是：

$$
\boldsymbol\chi_{p}^{A}(\mathbf r)
=
\mathcal N^{\mathrm{contr}}
\begin{bmatrix}
\mathbf b_1^{A,p}(\mathbf r) &
\mathbf b_2^{A,p}(\mathbf r) &
\cdots &
\mathbf b_M^{A,p}(\mathbf r)
\end{bmatrix}
\begin{bmatrix}
c_1 \\
c_2 \\
\vdots \\
c_M
\end{bmatrix}
$$

这里中间那个矩阵的形状是 $3\times M$：行对应 $p_x,p_y,p_z$ 三个角向分量（当前是以p壳层为例），列对应不同 primitive exponent。右边的 coefficient vector 形状是 $M\times 1$，乘出来就是 $3\times 1$ 的 contracted p shell 向量。

如果有多个 contraction，也就是 coefficient 有多列，那么只要把右边的系数列向量换成 coefficient matrix：

$$
\boldsymbol\chi_{p}^{A}(\mathbf r)
=
\begin{bmatrix}
\mathbf b_1^{A,p}(\mathbf r) &
\mathbf b_2^{A,p}(\mathbf r) &
\cdots &
\mathbf b_M^{A,p}(\mathbf r)
\end{bmatrix}
\mathbf C
$$

其中

$$
\boldsymbol\chi_{p}^{A}(\mathbf r)
\in \mathbb R^{3\times n_{\mathrm{ctr}}},
\qquad
\mathbf C\in \mathbb R^{M\times n_{\mathrm{ctr}}}
$$

结果矩阵的每一列是一组 contracted p shell：

$$
\boldsymbol\chi_{p}^{A}[:,k]
=
\begin{bmatrix}
\chi_{k,p_x}^{A} \\
\chi_{k,p_y}^{A} \\
\chi_{k,p_z}^{A}
\end{bmatrix}
$$

一般化到任意角动量 $l$，可以记成：

$$
\boldsymbol\chi_l^A(\mathbf r)
= \mathcal B_l^A(\mathbf r)\mathbf C
$$

其中

$$
\mathcal B_l^A(\mathbf r)
\in \mathbb R^{d_l\times n_{\mathrm{prim}}},
\qquad
\boldsymbol\chi_l^A(\mathbf r)
\in \mathbb R^{d_l\times n_{\mathrm{ctr}}}
$$

$d_l$ 是这个 shell 的角向函数数量：spherical 下 $d_l=2l+1$，Cartesian 下 $d_l=(l+1)(l+2)/2$。因此矩阵语言里的核心关系就是：

```text
primitive shell matrix  ×  contraction coefficient matrix  =  contracted shell matrix
(d_l × nprim)           ×  (nprim × nctr)                  =  (d_l × nctr)
```

再举一个对比：如果是 d shell，`l=2` 且 `nctr=1`，它也只有一个 contracted radial function，但 spherical 下会产生 5 个 d 型 AO，Cartesian 下会产生 6 个 d 型 AO。这里的“1”来自 coefficient 列数；“5 或 6”来自角向部分。

如果一个 shell 有多列 coefficient，那么多个 contraction 共用同一组 exponent。矩阵写法是：

$$
\mathbf B_{lm}^{A}(\mathbf r)
=
\begin{bmatrix}
B_{1lm}^{A}(\mathbf r) & B_{2lm}^{A}(\mathbf r) & \cdots & B_{n_{\mathrm{prim}}lm}^{A}(\mathbf r)
\end{bmatrix}
$$

其中每个 $B_{plm}^{A}$ 是一个已经带有角向分量的 primitive basis function：

$$
B_{plm}^{A}(\mathbf r)
= N_l(\alpha_p)r_A^le^{-\alpha_pr_A^2}Y_{lm}(\hat{\mathbf r}_A)
$$

于是对固定的 $l,m$，多个 contracted AO 可以写成：

$$
\boldsymbol\chi_{lm}^{A}(\mathbf r)
=
\begin{bmatrix}
\chi_{1lm}^{A}(\mathbf r) & \chi_{2lm}^{A}(\mathbf r) & \cdots & \chi_{n_{\mathrm{ctr}}lm}^{A}(\mathbf r)
\end{bmatrix}
= \mathbf B_{lm}^{A}(\mathbf r)\,\mathbf C
$$

这里 $\mathbf B_{lm}^{A}$ 是一个 $1\times n_{\mathrm{prim}}$ 的“primitive 函数行向量”，$\mathbf C$ 是一个 $n_{\mathrm{prim}}\times n_{\mathrm{ctr}}$ 的 coefficient matrix，乘出来就是 $1\times n_{\mathrm{ctr}}$ 的 contracted function 行向量。

$$
\mathbf C =
\begin{bmatrix}
 c_{11} & c_{12} & \cdots & c_{1n_{\mathrm{ctr}}} \\
 c_{21} & c_{22} & \cdots & c_{2n_{\mathrm{ctr}}} \\
 \vdots & \vdots & \ddots & \vdots \\
 c_{n_{\mathrm{prim}}1} & c_{n_{\mathrm{prim}}2} & \cdots & c_{n_{\mathrm{prim}}n_{\mathrm{ctr}}}
\end{bmatrix}
$$

如果暂时只关心径向部分，也可以把 $\mathbf B$ 理解成 primitive radial function 的行向量：

$$
\mathbf B_{l}^{A}(r_A)
=
\begin{bmatrix}
N_l(\alpha_1)r_A^le^{-\alpha_1r_A^2} &
N_l(\alpha_2)r_A^le^{-\alpha_2r_A^2} &
\cdots
\end{bmatrix}
$$

这时 $\mathbf B_l^A\mathbf C$ 得到的是多个 contracted radial functions；随后每个 radial function 再乘上同一个 shell 的各个角向函数，才得到最终 AO。

PySCF 代码中 `cs.shape == (nprim, nctr)`，正是这个 coefficient matrix 的形状。随后 `make_bas_env` 会把 primitive normalization 吸收到 `cs` 中，并把 exponent vector 与 coefficient matrix 存进 `_env`。

## 5. PySCF 内部 basis list 的数字含义

PySCF 的 orbital basis 内部格式大致是：

```python
[
    [l,
        [alpha_1, c_11, c_12, ...],
        [alpha_2, c_21, c_22, ...],
        ...
    ],
    ...
]
```

逐项对应：

| 位置 | 例子 | 含义 | 进入公式的位置 |
| --- | --- | --- | --- |
| `l` | `0` | angular momentum，`0=s`, `1=p`, `2=d` | $r^lY_{lm}$ 或 $x^ay^bz^c$ 中的总角动量 |
| `alpha_p` | `3.42525091` | 第 $p$ 个 primitive 的 Gaussian exponent | $e^{-\alpha_p r^2}$ |
| `c_pk` | `0.15432897` | 第 $p$ 个 primitive 对第 $k$ 个 contraction 的系数 | $\sum_p c_{pk}g_p$ |
| 行数 | 3 行 | primitive 数量 `nprim` | $p=1,\cdots,n_{\mathrm{prim}}$ |
| coefficient 列数 | 1 列或多列 | contraction 数量 `nctr` | $k=1,\cdots,n_{\mathrm{ctr}}$ |

非相对论计算中常见的是 `[l, [alpha, coeff...], ...]`。有时你会看到 `[l, kappa, [alpha, coeff...], ...]`，其中 `kappa` 用于 spinor/相对论相关的 basis 分支。第一周读 `mole.py` 时可以先把 `kappa=0` 当作普通情况。

## 6. H 的 STO-3G 例子

PySCF 自带的 `sto-3g.dat` 里，H 的一段是：

```text
H    S
      3.42525091             0.15432897
      0.62391373             0.53532814
      0.16885540             0.44463454
```

被解析成 PySCF 内部格式后，可以理解为：

```python
[[0,
  [3.42525091, 0.15432897],
  [0.62391373, 0.53532814],
  [0.16885540, 0.44463454]]]
```

这表示：

- `0`：这是一个 s shell，$l=0$。
- 三行数字：有 3 个 primitive。
- 每行第一列：$\alpha_1,\alpha_2,\alpha_3$。
- 每行第二列：$c_{11},c_{21},c_{31}$。
- 只有 1 列 coefficient，所以 `nctr = 1`。

公式就是：

$$
\chi_{1s}^{\mathrm H}(\mathbf r)
= \mathcal N^{\mathrm{contr}}
\sum_{p=1}^{3} c_p\,N_0(\alpha_p)e^{-\alpha_p r_H^2}
$$

代入数字写成：

$$
\chi_{1s}^{\mathrm H}(\mathbf r)
= \mathcal N^{\mathrm{contr}}\left[
0.15432897\,N_0(3.42525091)e^{-3.42525091r_H^2}
+0.53532814\,N_0(0.62391373)e^{-0.62391373r_H^2}
+0.44463454\,N_0(0.16885540)e^{-0.16885540r_H^2}
\right]
$$

这里 STO-3G 的意思可以直观看成：用 3 个 Gaussian primitive 去拟合一个 Slater-like 的 1s 原子轨道形状。

## 7. `SP` 行是什么意思

在 NWChem/Gaussian 风格基组里，有时会看到 `SP`，例如 STO-3G 中 C 的价层：

```text
C    SP
      2.9412494             -0.09996723             0.15591627
      0.6834831              0.39951283             0.60768372
      0.2222899              0.70011547             0.39195739
```

这不是一种新的角动量。它是一种紧凑写法，表示同一组 exponent 同时生成一个 s shell 和一个 p shell：

```python
# s shell
[0,
  [2.9412494, -0.09996723],
  [0.6834831,  0.39951283],
  [0.2222899,  0.70011547]]

# p shell
[1,
  [2.9412494, 0.15591627],
  [0.6834831, 0.60768372],
  [0.2222899, 0.39195739]]
```

第一列仍然是 exponent。第二列给 s contraction coefficient，第三列给 p contraction coefficient。

p shell 的 radial contraction 只有一组，但因为 $l=1$，球谐形式下会产生三个 AO：

$$
\chi_{p,m}^{A}(\mathbf r)
= \mathcal N^{\mathrm{contr}}
\sum_p c_pN_1(\alpha_p)r_Ae^{-\alpha_pr_A^2}Y_{1m}(\hat{\mathbf r}_A),
\qquad m=-1,0,1
$$

如果看 Cartesian 语言，就是同一个径向部分分别乘以 $x_A,y_A,z_A$。

## 8. 多列 coefficient 的例子

PySCF 的 parser 可以把共用 exponent 的多个 contraction 合并成 general contraction。例子：

```python
[[0,
  [13.6267, 0.17523, 0.0],
  [1.99935, 0.893483, 0.0],
  [0.382993, 0.0, 1.0]]]
```

这表示一个 s shell，`nprim = 3`，`nctr = 2`。两个 contracted s 函数分别是：

$$
\chi_1(\mathbf r)
= \sum_{p=1}^3 c_{p1}N_0(\alpha_p)e^{-\alpha_pr^2}
$$

$$
\chi_2(\mathbf r)
= \sum_{p=1}^3 c_{p2}N_0(\alpha_p)e^{-\alpha_pr^2}
$$

矩阵形式是：

$$
\begin{bmatrix}
\chi_1 & \chi_2
\end{bmatrix}
=
\begin{bmatrix}
B_1 & B_2 & B_3
\end{bmatrix}
\begin{bmatrix}
0.17523 & 0.0 \\
0.893483 & 0.0 \\
0.0 & 1.0
\end{bmatrix}
$$

其中

$$
B_p(\mathbf r)=N_0(\alpha_p)e^{-\alpha_pr^2}
$$

## 9. 从 `basis` 到 `_bas` 和 `_env`

`Mole.build()` 读入 basis 后，大致会走这条路：

```text
mol.basis
  -> _parse_default_basis
  -> format_basis
  -> make_env
  -> make_bas_env
  -> mol._bas, mol._env
```

其中 `format_basis` 负责把字符串、列表、混合写法转换成统一内部格式：

```python
{atom: [[l, [alpha, c1, c2, ...], ...], ...]}
```

之后 `make_bas_env` 会对每个 shell 做几件事：

```python
angl = b[0]
es = b_coeff[:, 0]
cs = b_coeff[:, 1:]
nprim, nctr = cs.shape
cs = cs * gto_norm(angl, es)
```

也就是说：

- `es` 是 exponent vector，形状是 $(n_{\mathrm{prim}},)$。
- `cs` 是 coefficient matrix，形状是 $(n_{\mathrm{prim}}, n_{\mathrm{ctr}})$。
- PySCF 会把 primitive normalization 吸收到 coefficient matrix 里。
- 如果 `NORMALIZE_GTO = True`，还会继续归一化 contracted AO。

最后 `_bas` 的每一行是一个 shell 的索引信息：

```text
[atom_id, l, nprim, nctr, kappa, ptr_exp, ptr_coeff, 0]
```

对应表：

| `_bas` 列 | 含义 |
| --- | --- |
| `atom_id` | 这个 shell 属于哪个原子 |
| `l` | 角动量 |
| `nprim` | primitive 数量 |
| `nctr` | contraction 数量 |
| `kappa` | spinor/相对论相关参数，普通非相对论通常为 0 |
| `ptr_exp` | exponent 在 `_env` 里的起始位置 |
| `ptr_coeff` | coefficient 在 `_env` 里的起始位置 |
| 最后一列 | 保留位 |

`_env` 则是真正存浮点数的长数组。对每个 shell，它先存 exponent，再存归一化后的 coefficient。

所以可以这样理解：

```text
basis list 适合人读
_bas       适合积分库定位每个 shell 的元信息
_env       适合积分库存取 exponent 和 coefficient 的连续内存
```

## 10. shell 数和 AO 数不要混淆

一个 shell 不是一个 AO。一个 shell 会展开成多少 AO，取决于 $l$、`nctr` 和 spherical/cartesian 选择。

球谐 AO 数：

$$
N_{\mathrm{AO}}^{\mathrm{sph}}(\text{shell}) = (2l+1)n_{\mathrm{ctr}}
$$

Cartesian AO 数：

$$
N_{\mathrm{AO}}^{\mathrm{cart}}(\text{shell}) = \frac{(l+1)(l+2)}{2}n_{\mathrm{ctr}}
$$

例如一个 p shell，`l=1, nctr=1`，不是 1 个 AO，而是 3 个 AO。一个 d shell，`l=2, nctr=1`，球谐下是 5 个 AO，Cartesian 下是 6 个 AO。

## 11. 读源码时建议盯住哪里

先按这个顺序读，不用一次读完所有基组解析器：

1. `pyscf/gto/mole.py` 的 `format_basis`：看 PySCF 统一内部 basis list 的目标格式。
2. `pyscf/gto/basis/parse_nwchem.py` 的 `_parse`：看文本行如何变成 `[l, [alpha, coeff...]]`。
3. `pyscf/gto/mole.py` 的 `gto_norm`：看 primitive normalization 的公式。
4. `pyscf/gto/mole.py` 的 `make_bas_env`：看 `alpha` 和 `coeff` 如何进入 `_bas/_env`。
5. `pyscf/gto/moleintor.py` 的 `make_loc`：看一个 shell 如何对应到 AO 数量。

## 12. MYNOTE: 我现在应该掌握到什么程度

第一周不需要把所有 basis set 家族背下来。更重要的是看到一行数字时能说出：

```text
第一列是 exponent alpha，控制 e^{-alpha r^2}
后面的列是 contraction coefficients c_pk
外层的 l 决定 s/p/d/f 和角向函数数量
行数决定 nprim
coefficient 列数决定 nctr
```

再进一步，要能把这个列表：

```python
[1,
 [2.0, 0.3],
 [0.5, 0.7]]
```

翻译成一句话：这是一个 p shell，两个 primitive，一个 contraction；它先形成一个径向收缩函数，再乘以三个 p 方向的角向函数，因此 spherical/cartesian 下都会贡献 3 个 AO。
