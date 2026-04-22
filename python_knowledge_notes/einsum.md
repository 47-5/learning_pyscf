# Python 知识点：`numpy.einsum`

## 一句话理解

`numpy.einsum` 可以把“按下标求和”的规则直接写在字符串里。

最核心的思想只有两条：

1. 同一个下标如果在输入里重复出现而没有出现在输出里，就对它求和。
2. 输出里下标的顺序，决定结果数组各个轴的顺序。

## 本文默认约定

为了和你更习惯的线性代数语言一致，本文默认：

- 单个向量默认写成列向量。
- 单个点的坐标写成

$$
\mathbf r =
\begin{bmatrix}
x \\
y \\
z
\end{bmatrix}
$$

但要注意，NumPy 代码里常常会把很多个坐标按“行”堆起来存。

例如 `format_atom` 里的：

```python
c = numpy.array([a[1] for a in fmt_atoms], dtype=numpy.double)
```

如果有 $N$ 个原子，那么 `c.shape == (N, 3)`，也就是：

```text
第 0 维: 第几个原子
第 1 维: x/y/z 三个坐标分量
```

所以：

- 数学上看单个原子时，我会写列向量。
- 看源码里的 `c` 时，要记得它是“很多个坐标行向量堆起来的矩阵”。

## 一个最简单的例子：点积

```python
numpy.einsum('i,i->', a, b)
```

表示：

$$
\sum_i a_i b_i
$$

因为：

- `a` 的下标是 `i`
- `b` 的下标也是 `i`
- 输出里没有 `i`

所以 `i` 被求和掉，结果是一个标量。

## 矩阵-向量乘法

```python
numpy.einsum('ij,j->i', A, x)
```

表示：

$$
y_i = \sum_j A_{ij}x_j
$$

也就是普通的矩阵乘法：

$$
\mathbf y = A\mathbf x
$$

这里：

- `A` 的形状是 `(i, j)`
- `x` 的形状是 `(j,)`
- 输出是 `(i,)`

## 矩阵-矩阵乘法

```python
numpy.einsum('ij,jk->ik', A, B)
```

表示：

$$
C_{ik} = \sum_j A_{ij}B_{jk}
$$

也就是：

$$
C = AB
$$

## 看 `format_atom` 里的这一句

`pyscf/gto/mole.py` 里有：

```python
c = numpy.einsum('ix,kx->ki', axes * unit, c - origin)
```

先看两个输入张量：

```python
axes * unit
```

形状是 `(3, 3)`，标记成 `ix`。

```python
c - origin
```

如果有 $N$ 个原子，形状是 `(N, 3)`，标记成 `kx`。

这里各个下标的含义可以这样读：

- `k`：第几个原子
- `x`：旧坐标系中的分量索引，取值为 $x,y,z$
- `i`：新坐标系中的分量索引，取值为 $x',y',z'$

输出是 `ki`，所以结果形状是 `(N, 3)`。

而中间重复出现、但没有出现在输出里的 `x`，说明要对旧坐标分量求和。

因此它等价于：

$$
c_{ki}^{\mathrm{new}} = \sum_x (\mathrm{axes}_{ix}\cdot \mathrm{unit})(c_{kx}-\mathrm{origin}_x)
$$

## 用列向量语言怎么写

先只看第 $k$ 个原子。

记它原来的坐标列向量是：

$$
\mathbf r_k =
\begin{bmatrix}
x_k \\
y_k \\
z_k
\end{bmatrix}
$$

原点平移向量是：

$$
\mathbf r_0 =
\begin{bmatrix}
o_x \\
o_y \\
o_z
\end{bmatrix}
$$

坐标轴变换矩阵记为：

$$
A = \mathrm{axes}
$$

那么这句 `einsum` 对单个原子的数学含义就是：

$$
\mathbf r_k' = (\mathrm{unit}\cdot A)(\mathbf r_k - \mathbf r_0)
$$

也就是：

1. 先减去新原点
2. 再投影到新坐标轴上
3. 再把单位统一换成 Bohr

## 为什么代码里不是直接写 `A @ r`

因为代码里的 `c` 不是一个列向量，而是很多原子的坐标按行排成的矩阵：

$$
C =
\begin{bmatrix}
- & \mathbf r_1^T & - \\
- & \mathbf r_2^T & - \\
& \vdots & \\
- & \mathbf r_N^T & -
\end{bmatrix}
\in \mathbb R^{N\times 3}
$$

如果按列向量习惯，把所有原子写成一列一列排开，则更自然的总公式是：

$$
R =
\begin{bmatrix}
\mathbf r_1 & \mathbf r_2 & \cdots & \mathbf r_N
\end{bmatrix}
\in \mathbb R^{3\times N}
$$

那么整体变换写成：

$$
R' = (\mathrm{unit}\cdot A)(R - \mathbf r_0\mathbf 1^T)
$$

这里 $\mathbf 1$ 是长度为 $N$ 的全 1 列向量。

但是 NumPy 代码里用的是按行存储的 `c`，所以等价地可以写成：

$$
C' = (C - \mathbf 1\mathbf r_0^T)(\mathrm{unit}\cdot A)^T
$$

这正是 `einsum('ix,kx->ki', ...)` 在做的事。

## 它和普通矩阵乘法的等价形式

这一句：

```python
c = numpy.einsum('ix,kx->ki', axes * unit, c - origin)
```

完全等价于：

```python
c = (c - origin) @ (axes * unit).T
```

或者如果只看单个原子的列向量语言：

```text
r_new = (axes * unit) @ (r_old - origin)
```

两种写法没有本质区别，只是“很多点按行存”与“单个点按列写”的视角不同。

## 怎么快速读 `einsum`

以后看到 `einsum`，你可以按这个顺序读：

1. 先看每个输入的下标各代表什么物理量。
2. 找出重复下标里哪些没有出现在输出里。
3. 这些下标就是要求和掉的维度。
4. 输出下标的顺序，就是结果数组轴的顺序。

例如：

```python
'ix,kx->ki'
```

就可以读成：

```text
对旧坐标分量 x 求和
输入: 坐标轴矩阵 (i,x) 和原子坐标 (k,x)
输出: 每个原子在新坐标系下的坐标 (k,i)
```

## MYNOTE：现在至少要掌握什么

读 PySCF 源码时，不需要把 `einsum` 当成神秘黑魔法。对这一行来说，你只要能把它翻译成：

```text
先平移，再做坐标轴变换，再换单位
```

并且知道它等价于矩阵乘法，就已经足够了。