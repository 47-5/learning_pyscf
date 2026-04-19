# Symmetry-adapted AO basis 概念笔记

## 一句话理解

原始 AO 基函数通常以各个原子为中心定义。它们很适合做积分和表达局域化学直觉，但单个 AO 不一定具有整个分子点群下的明确对称性。

启用 symmetry 后，PySCF 可以把这些原始 AO 做线性组合，得到一组按不可约表示（irreducible representation, irrep）分类的 symmetry-adapted AO basis。

重要的是：

> symmetry-adapted AO basis 没有改变原始 basis set 所张成的计算空间，只是在同一个空间里换了一套更适合利用分子对称性的坐标。

## 关于“空间”的理解

这里的“空间”不是三维物理空间，而是由基函数张成的函数向量空间。

假设原始 AO basis 是：

$$
\{\chi_1(\mathbf r), \chi_2(\mathbf r), \ldots, \chi_N(\mathbf r)\}
$$

它们的所有线性组合构成一个 AO 空间：

$$
V_{\mathrm{AO}} = \mathrm{span}\{\chi_1, \chi_2, \ldots, \chi_N\}.
$$

任何一个分子轨道都可以写成这个空间里的一个向量：

$$
\psi(\mathbf r) = \sum_{\mu=1}^N c_\mu \chi_\mu(\mathbf r).
$$

矩阵形式可以把所有 AO 排成一个“行向量”：

$$
\boldsymbol{\chi}(\mathbf r)
= \begin{bmatrix}
\chi_1(\mathbf r) & \chi_2(\mathbf r) & \cdots & \chi_N(\mathbf r)
\end{bmatrix},
$$

于是分子轨道可以写成：

$$
\psi(\mathbf r) = \boldsymbol{\chi}(\mathbf r)\mathbf c,
$$

其中：

$$
\mathbf c = \begin{bmatrix}
c_1 \\
c_2 \\
\vdots \\
c_N
\end{bmatrix}.
$$

启用 symmetry 后，程序构造新的基函数：

$$
\{\phi_1(\mathbf r), \phi_2(\mathbf r), \ldots, \phi_N(\mathbf r)\}.
$$

每个 $\phi_i$ 都是原始 AO 的线性组合：

$$
\phi_i(\mathbf r) = \sum_{\mu=1}^N \chi_\mu(\mathbf r) C_{\mu i}.
$$

矩阵形式为：

$$
\boldsymbol{\phi}(\mathbf r) = \boldsymbol{\chi}(\mathbf r) C,
$$

其中 $C$ 是从原始 AO 到新基函数的变换矩阵。如果这组 $\phi_i$ 仍然张成同一个空间：

$$
\mathrm{span}\{\phi_1, \ldots, \phi_N\}
=
\mathrm{span}\{\chi_1, \ldots, \chi_N\},
$$

那么它不是换了 basis set，而只是换了这个向量空间的坐标系。

类比线性代数：二维平面可以用 $x/y$ 坐标表示，也可以用旋转 $45^\circ$ 后的坐标表示。平面没有变，只是坐标轴变了。

在这里：

| 对象 | 类比 |
|---|---|
| 原始 AO basis | 一套以原子为中心的坐标轴 |
| symmetry-adapted AO basis | 一套按分子对称性整理后的坐标轴 |
| AO function space | 两者共同张成的计算空间 |


## NAO 是什么

在 PySCF 的 `gto/mole.py` 语境里，`nao` 通常表示：

> number of atomic orbitals，也就是 AO 基函数的数量。

如果原始 AO basis 是：

$$
\{\chi_1, \chi_2, \ldots, \chi_N\},
$$

那么：

$$
N = N_{\mathrm{AO}} = \mathrm{nao} = \dim V_{\mathrm{AO}}.
$$

也就是说，`nao` 是 AO 函数空间的维度。

在 PySCF 中常见写法包括：

```python
mol.nao
mol.nao_nr()
nao
```

其中 `nao_nr()` 里的 `nr` 是 non-relativistic，即非相对论 AO 数。普通分子 RHF/RKS 计算里，通常可以先把 `mol.nao` 和 `mol.nao_nr()` 近似看成同一个概念。

如果：

```python
nao = mol.nao_nr()
```

那么常见矩阵和张量的形状是：

| 对象 | 数学含义 | 常见 shape |
|---|---|---|
| overlap matrix | $S_{\mu\nu}=\langle\chi_\mu|\chi_\nu\rangle$ | `(nao, nao)` |
| core Hamiltonian | $h_{\mu\nu}$ | `(nao, nao)` |
| density matrix | $D_{\mu\nu}$ | `(nao, nao)` |
| Fock matrix | $F_{\mu\nu}$ | `(nao, nao)` |
| MO coefficients | $C_{\mu i}$ | `(nao, nmo)` |
| two-electron integrals | $(\mu\nu|\lambda\sigma)$ | `(nao, nao, nao, nao)` in full form |

例如 H2O/STO-3G 中，可以粗略理解为：

| 原子 | AO 数 |
|---|---:|
| O | 5, namely $1s,2s,2p_x,2p_y,2p_z$ |
| H | 1, namely $1s$ |
| H | 1, namely $1s$ |
| total | 7 |

所以 H2O/STO-3G 的 `nao` 通常是 7。

注意：在某些量子化学文献中，NAO 也可能表示 natural atomic orbital。但在 `mole.py` 这类 PySCF 源码上下文里，`nao` 基本表示 number of atomic orbitals。

可以写进源码注释：

```python
# MYNOTE: nao = number of AO basis functions, i.e. dim(V_AO).
# MYNOTE: one-electron AO matrices usually have shape (nao, nao).
# MYNOTE: full ERI tensor has shape (nao, nao, nao, nao).
```

## irreps 是什么

每个点群都有自己的一套不可约表示（irreducible representations, irreps）。例如：

| 点群 | 常见 irreps |
|---|---|
| $C_{2v}$ | $A_1, A_2, B_1, B_2$ |
| $D_{2h}$ | $A_g, B_{1g}, B_{2g}, B_{3g}, A_u, B_{1u}, B_{2u}, B_{3u}$ |
| $C_s$ | $A', A''$ |
| $C_i$ | $A_g, A_u$ |

irrep 可以先直观理解为：

> 一个轨道或基函数组合在点群操作下的“对称性类型”。

比如 $C_{2v}$ 中的 $A_1$ 通常代表完全对称的一类函数；$B_1$、$B_2$ 则代表在某些旋转或反射操作下会变号的函数类型。严格含义要看该点群的 character table。

从线性代数角度看，点群操作在 AO 空间中形成一个 representation。如果通过某个变换矩阵 $C$，所有对称操作矩阵都能被同时变成**块**对角形式：

$$
C^{-1}D(g)C
=
\begin{bmatrix}
D^{(\Gamma_1)}(g) & 0 & \cdots \\
0 & D^{(\Gamma_2)}(g) & \cdots \\
\vdots & \vdots & \ddots
\end{bmatrix},
\quad g\in G,
$$

那么每个不可再分的小块 $D^{(\Gamma)}(g)$ 就对应一个 irrep $\Gamma$。原始 AO 空间可以被分解成不同 irrep 的子空间：

$$
V_{\mathrm{AO}} = \bigoplus_\Gamma V_\Gamma.
$$

这里的 $\oplus$ 表示直和，也就是不同对称性子空间之间不会互相混合。

注意：某个点群有一套 irreps，但某个具体分子和具体 basis set 不一定会出现所有 irreps。例如一个 $C_{2v}$ 分子理论上有 $A_1,A_2,B_1,B_2$，但某个很小的 basis 里可能没有 $A_2$ 类型的 symmetry-adapted AO。

## 点群操作作用在 AO basis 空间上是什么意思

点群操作包括旋转、反射、反演等。它们会改变空间坐标，例如：

$$
\mathbf r \mapsto g\mathbf r.
$$

当这个操作作用在函数上时，可以写成：

$$
(\hat R_g f)(\mathbf r) = f(g^{-1}\mathbf r).
$$

这里 $\hat R_g$ 是点群操作 $g$ 对函数的作用算符。对 AO 基函数而言：

$$
\hat R_g \chi_\mu(\mathbf r)
= \sum_{\nu=1}^N \chi_\nu(\mathbf r)D_{\nu\mu}(g).
$$

矩阵形式为：

$$
\hat R_g\boldsymbol{\chi}(\mathbf r)
= \boldsymbol{\chi}(\mathbf r)D(g).
$$

其中：

| 符号 | 含义 |
|---|---|
| $g$ | 一个点群操作，例如反射、旋转、反演 |
| $G$ | 整个点群，也就是所有对称操作的集合 |
| $\hat R_g$ | 点群操作 $g$ 对函数的作用算符 |
| $D(g)$ | 操作 $g$ 在 AO basis 空间中的矩阵表示 |

这就是 representation 的来源：点群操作可以被表示成作用在 AO 空间上的矩阵。

直观地说：

- 如果对称操作把一个原子映射到另一个等价原子，这个原子上的 AO 会被映射到另一个原子上的 AO。
- 如果原子中心不动，但轨道有方向性，例如 $p_x,p_y,p_z$，那么对称操作可能让它变号，或者变成另一个方向的 p 轨道。
- 对于 s 轨道，方向性较弱，很多时候主要表现为在等价原子之间互相交换。

## 最简单例子：两个等价 H 的 1s

考虑 H$_2$或任意有两个等价 H 的片段。原始 AO 是：

$$
\chi_A = \text{左边 H 的 }1s,
\qquad
\chi_B = \text{右边 H 的 }1s.
$$

某个对称操作 $s$ 会交换两个 H：

$$
\hat R_s\chi_A = \chi_B,
\qquad
\hat R_s\chi_B = \chi_A.
$$

在原始 AO basis $(\chi_A,\chi_B)$ 中，这个操作的矩阵可以写成：

$$
D(s) =
\begin{bmatrix}
0 & 1 \\
1 & 0
\end{bmatrix}.
$$

恒等操作的矩阵为：

$$
D(E) =
\begin{bmatrix}
1 & 0 \\
0 & 1
\end{bmatrix}.
$$

如果改用线性组合：

$$
\phi_+ = \frac{1}{\sqrt 2}(\chi_A + \chi_B),
\qquad
\phi_- = \frac{1}{\sqrt 2}(\chi_A - \chi_B),
$$

再做交换操作：

$$
\hat R_s\phi_+ = +\phi_+,
\qquad
\hat R_s\phi_- = -\phi_-.
$$

这说明 $\phi_+$ 和 $\phi_-$ 已经是对称性很明确的组合。它们分别属于不同的 irreps。

把这个变换写成矩阵：

$$
C_{\mathrm{sym}} = \frac{1}{\sqrt 2}
\begin{bmatrix}
1 & 1 \\
1 & -1
\end{bmatrix},
\qquad
\begin{bmatrix}
\phi_+ & \phi_-
\end{bmatrix}
=
\begin{bmatrix}
\chi_A & \chi_B
\end{bmatrix}C_{\mathrm{sym}}.
$$

在新 basis 下，交换操作变成对角矩阵：

$$
C_{\mathrm{sym}}^{-1}D(s)C_{\mathrm{sym}}
=
\begin{bmatrix}
1 & 0 \\
0 & -1
\end{bmatrix}.
$$

这个矩阵版本很重要：它说明 symmetry-adapted basis 的作用就是把“会互相交换/混合”的原始 AO 重新组合成“对称操作下只乘以简单因子”的基函数。

## 为什么原始 AO 不一定已经按对称性分好

原始 AO 是局域在原子上的，例如：

$$
\text{H}_{\mathrm{left}}1s,
\quad
\text{H}_{\mathrm{right}}1s,
\quad
\text{O }2s,
\quad
\text{O }2p_x,
\quad
\text{O }2p_y,
\quad
\text{O }2p_z.
$$

其中有些 AO 单独看可能已经有明确对称性，例如位于对称轴上的某些 $p$ 轨道。

但很多原始 AO，尤其是位于等价原子上的 AO，会在对称操作下互相交换。单独一个 $\text{H}_{\mathrm{left}}1s$ 并不是整个分子点群下的“纯净对称性类型”。

通过线性组合：

$$
\text{H}_{\mathrm{left}}1s + \text{H}_{\mathrm{right}}1s,
\qquad
\text{H}_{\mathrm{left}}1s - \text{H}_{\mathrm{right}}1s,
$$

就能得到明确属于某些 irreps 的组合。

## symmetry-adapted AO basis 是怎么构造的

概念流程如下：

1. 检测分子的点群，例如 $C_{2v}$、$D_{2h}$、$C_{\infty v}$。
2. 列出这个点群的对称操作 $g\in G$。
3. 判断每个对称操作如何作用在原始 AO basis 上。
4. 得到每个操作在 AO 空间中的矩阵表示 $D(g)$。
5. 使用点群的 character table 和投影算符，把 AO 空间投影到各个 irrep 对应的子空间。
6. 对得到的线性组合归一化，并处理线性相关。
7. 得到按 irrep 分组的 symmetry-adapted AO basis。

对某个 irrep $\Gamma$，投影算符的函数形式可以写成：

$$
\hat P^{(\Gamma)}
=
\frac{d_\Gamma}{|G|}
\sum_{g\in G}
\chi^{(\Gamma)}(g)^*\hat R_g.
$$

矩阵形式为：

$$
P^{(\Gamma)}
=
\frac{d_\Gamma}{|G|}
\sum_{g\in G}
\chi^{(\Gamma)}(g)^*D(g).
$$

其中：

| 符号 | 含义 |
|---|---|
| $\Gamma$ | 某个 irrep，例如 $A_1$ 或 $B_{2u}$ |
| $d_\Gamma$ | 这个 irrep 的维数 |
| $|G|$ | 点群中对称操作的数量 |
| $\chi^{(\Gamma)}(g)$ | irrep $\Gamma$ 对操作 $g$ 的 character |
| $D(g)$ | 操作 $g$ 在原始 AO 空间中的矩阵表示 |
| $P^{(\Gamma)}$ | 把 AO 空间投影到 $\Gamma$ 子空间的投影矩阵 |

这里有两个不同的 $\chi$，容易混：$\chi_\mu(\mathbf r)$ 是 AO 基函数，$\chi^{(\Gamma)}(g)$ 是 character table 里的 character。上下文不同。

投影矩阵的作用可以写成：

$$
\mathbf v^{(\Gamma)} = P^{(\Gamma)}\mathbf v,
$$

其中 $\mathbf v$ 是原始 AO 空间中的某个系数向量，$\mathbf v^{(\Gamma)}$ 是它在 irrep $\Gamma$ 子空间中的部分。

## 和 PySCF 字段的对应

在 `pyscf/gto/mole.py` 中，PySCF 会调用：

```python
self.symm_orb, self.irrep_id = \
    symm.symm_adapted_basis(self, groupname, orig, axes)
```

然后：

```python
self.irrep_name = [symm.irrep_id2name(groupname, ir)
                   for ir in self.irrep_id]
```

可以这样理解：

| 字段 | 含义 |
|---|---|
| `mol.symm_orb` | 按 irrep 分组的 symmetry-adapted AO 系数矩阵 |
| `mol.symm_orb[i]` | 第 `i` 个 irrep 对应的 symmetry-adapted AO basis |
| `mol.irrep_name[i]` | 第 `i` 组 symmetry-adapted AO 的人类可读标签，例如 `A1`, `B1`, `Ag` |
| `mol.irrep_id[i]` | 第 `i` 组 symmetry-adapted AO 的 PySCF 内部整数编号 |

如果：

```python
mol.symm_orb[i].shape == (nao, n_i)
```

则表示这个 irrep 下有 $n_i$ 个 symmetry-adapted AO。

记：

$$
C_\Gamma = \texttt{mol.symm\_orb[i]},
\qquad
C_\Gamma \in \mathbb R^{N_{\mathrm{AO}}\times n_\Gamma}.
$$

那么这个 irrep 下的 symmetry-adapted AO 可以写成：

$$
\boldsymbol{\phi}_\Gamma(\mathbf r)
=
\boldsymbol{\chi}(\mathbf r)C_\Gamma.
$$

如果把所有 irrep 的系数矩阵横向拼起来：

$$
C_{\mathrm{sym}}
=
\begin{bmatrix}
C_{\Gamma_1} & C_{\Gamma_2} & \cdots & C_{\Gamma_m}
\end{bmatrix},
$$

那么完整 symmetry-adapted AO basis 为：

$$
\boldsymbol{\phi}(\mathbf r)
=
\boldsymbol{\chi}(\mathbf r)C_{\mathrm{sym}}.
$$

## 矩阵为什么会按 irrep 分块

如果某个矩阵对应的算符和所有点群操作都对易，例如 Fock 矩阵对应的有效哈密顿量保持分子点群对称性，那么不同 irreps 之间不会互相耦合。

用 symmetry-adapted AO basis 表示 Fock 矩阵，可以写成：

$$
F^{\mathrm{sym}}
=
C_{\mathrm{sym}}^\dagger F C_{\mathrm{sym}}.
$$

如果 basis 非正交，更严格的程序实现会同时考虑 AO overlap 矩阵 $S$。概念上，关键结论是不同 irreps 之间的 block 为零：

$$
F_{\Gamma\Lambda}
=
C_\Gamma^\dagger F C_\Lambda
=0,
\qquad
\Gamma\ne\Lambda.
$$

于是矩阵具有块对角结构：

$$
F^{\mathrm{sym}}
=
\begin{bmatrix}
F_{\Gamma_1} & 0 & \cdots \\
0 & F_{\Gamma_2} & \cdots \\
\vdots & \vdots & \ddots
\end{bmatrix}.
$$

这就是 symmetry 能让计算和轨道分析更清楚的根本原因：不同对称性类型的轨道不会随便混在一起。

## “不是换了基组，而是换了坐标”

这点非常重要。

输入的 basis set 仍然是实际计算的 basis set。比如：

```python
basis = 'def2-svp'
```

这个 basis set 定义了 AO 空间本身。启用 symmetry 后，PySCF 只是把这个 AO 空间重新表达成按 irreps 分类的线性组合。

所以不应该说：

> 原始基组不是真正用于计算的基组。

更准确的说法是：

> 原始 basis set 定义了计算空间；symmetry-adapted AO basis 是同一个空间里的另一套更适合利用点群对称性的基函数坐标。

## 可以写进源码注释的 MYNOTE

```python
# MYNOTE: The input AO basis defines the function space used in the calculation.
# MYNOTE: Symmetry does not change this AO space; it builds a new coordinate system inside it.
# MYNOTE: symm_orb stores linear combinations of original AOs that transform as irreps.
# MYNOTE: A point-group operation acts on AO space by exchanging, changing signs, or mixing AO functions.
# MYNOTE: In matrix language, symm_orb[i] is a coefficient matrix C_Gamma for one irrep subspace.
# MYNOTE: symmetry-adapted AOs allow matrices and orbitals to be organized by irreducible representations.
```

## 第一周掌握目标

目前只需要掌握：

1. 原始 AO basis 是以原子为中心定义的。
2. 点群操作会在 AO 空间中交换、变号或混合 AO。
3. 每个点群有自己的一套 irreps。
4. irrep 是基函数或轨道在点群操作下的对称性类型。
5. symmetry-adapted AO basis 是原始 AO 的线性组合。
6. 它们和原始 AO 张成同一个函数空间。
7. 这些线性组合可以按 irreps 分组。
8. PySCF 用 `mol.symm_orb` 保存这些线性组合。