# `format_atom` 里的平移、换轴与旋转

## 背景代码

`pyscf/gto/mole.py` 的 `format_atom` 末尾有：

```python
if axes is None:
    axes = numpy.eye(3)

unit = _length_in_au(unit)
c = numpy.array([a[1] for a in fmt_atoms], dtype=numpy.double)
c = numpy.einsum('ix,kx->ki', axes * unit, c - origin)
z = [a[0] for a in fmt_atoms]
return list(zip(z, c.tolist()))
```

所以这部分代码做的事情可以概括成：

```text
读入原子坐标
    -> 以 origin 为新原点平移
    -> 用 axes 指定的新坐标轴重写坐标
    -> 把单位统一换成 Bohr
```

## 本文默认约定：向量写成列向量

为了和常见线性代数写法一致，本文默认单个坐标向量写成列向量：

$$
\mathbf r =
\begin{bmatrix}
x \\
y \\
z
\end{bmatrix}
$$

但源码里的 `c` 是许多原子坐标按行堆起来的 NumPy 数组：

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

这两种写法只是存储习惯不同，不改变几何含义。

## 第一步：平移

设新原点是：

$$
\mathbf r_0 =
\begin{bmatrix}
o_x \\
o_y \\
o_z
\end{bmatrix}
$$

那么一个原子坐标先变成：

$$
\tilde{\mathbf r} = \mathbf r - \mathbf r_0
$$

这一步没有旋转，只是把坐标原点改掉。

例如，如果原来某个原子在

$$
\mathbf r=
\begin{bmatrix}
1 \\
2 \\
3
\end{bmatrix}
$$

而新原点是

$$
\mathbf r_0=
\begin{bmatrix}
1 \\
1 \\
1
\end{bmatrix}
$$

那么平移后坐标就是

$$
\tilde{\mathbf r}=
\begin{bmatrix}
0 \\
1 \\
2
\end{bmatrix}
$$

## 第二步：`axes` 到底表示什么

这里最关键的一点是：

**`axes` 描述的是“新坐标轴在旧坐标系里的方向”。**

也就是说，如果写

$$
A = \mathrm{axes}
$$

那么它的三行分别是：

$$
A =
\begin{bmatrix}
(\mathbf e'_x)^T \\
(\mathbf e'_y)^T \\
(\mathbf e'_z)^T
\end{bmatrix}
$$

其中：

- $\mathbf e'_x$ 是新 x 轴在旧坐标系里的单位向量
- $\mathbf e'_y$ 是新 y 轴在旧坐标系里的单位向量
- $\mathbf e'_z$ 是新 z 轴在旧坐标系里的单位向量

于是，一个点在新坐标系下的三个分量就是对这些新轴做投影：

$$
x' = (\mathbf e'_x)^T\tilde{\mathbf r}
$$

$$
y' = (\mathbf e'_y)^T\tilde{\mathbf r}
$$

$$
z' = (\mathbf e'_z)^T\tilde{\mathbf r}
$$

合在一起就是：

$$
\mathbf r' = A\tilde{\mathbf r} = A(\mathbf r-\mathbf r_0)
$$

再乘上单位换算系数 `unit`：

$$
\mathbf r' = (\mathrm{unit}\cdot A)(\mathbf r-\mathbf r_0)
$$

## 这到底是“分子在转”还是“坐标轴在转”

这里更自然的理解是：

**坐标轴在变，分子本身并没有先在物理空间里动。**

也就是说：

- 原子还在原来的空间位置
- 你换了一套新的参考坐标轴
- 所以同一个点的坐标数值变了

这叫“换坐标系”或者“被动旋转”。

如果你更习惯“把分子转动一个角度”的想法，那叫“主动旋转”。对于正交旋转矩阵，这两种观点互为逆变换，也就是互相差一个转置。

## 一个最简单的例子：不旋转

如果

$$
A = I =
\begin{bmatrix}
1&0&0\\
0&1&0\\
0&0&1
\end{bmatrix}
$$

那么：

$$
\mathbf r' = \mathbf r - \mathbf r_0
$$

也就是没有换轴，只有平移。

这也对应源码中的：

```python
if axes is None:
    axes = numpy.eye(3)
```

## 一个交换坐标轴的例子

设

$$
A=
\begin{bmatrix}
0&1&0\\
1&0&0\\
0&0&1
\end{bmatrix}
$$

它表示：

- 新 x 轴沿旧 y 方向
- 新 y 轴沿旧 x 方向
- 新 z 轴不变

若某点相对新原点的坐标是

$$
\tilde{\mathbf r}=
\begin{bmatrix}
1\\
2\\
3
\end{bmatrix}
$$

则

$$
\mathbf r' = A\tilde{\mathbf r}
=
\begin{bmatrix}
0&1&0\\
1&0&0\\
0&0&1
\end{bmatrix}
\begin{bmatrix}
1\\
2\\
3
\end{bmatrix}
=
\begin{bmatrix}
2\\
1\\
3
\end{bmatrix}
$$

所以同一个点在新坐标系中被记成 $(2,1,3)$。

## 一个真正的“绕 z 轴换轴”例子

如果新坐标轴相对旧坐标轴发生了旋转，那么可以写成一个正交矩阵。

例如设：

$$
A=
\begin{bmatrix}
\cos\theta & \sin\theta & 0\\
-\sin\theta & \cos\theta & 0\\
0&0&1
\end{bmatrix}
$$

则

$$
\mathbf r' = A\mathbf r
$$

也就是

$$
x' = x\cos\theta + y\sin\theta
$$

$$
y' = -x\sin\theta + y\cos\theta
$$

$$
z' = z
$$

这表示：你不是把点先转了，而是把参考坐标轴换成了一套新的方向。

## 对很多原子一起写的矩阵形式

如果按列向量习惯，把所有原子坐标写成：

$$
R =
\begin{bmatrix}
\mathbf r_1 & \mathbf r_2 & \cdots & \mathbf r_N
\end{bmatrix}
\in \mathbb R^{3\times N}
$$

那么整体变换就是：

$$
R' = (\mathrm{unit}\cdot A)(R - \mathbf r_0\mathbf 1^T)
$$

这里：

- $\mathbf r_0\mathbf 1^T$ 把原点向量复制到每一列
- 左乘 $A$ 表示对每个原子都做同一个坐标轴变换

## 为什么源码里是 `N×3` 而不是 `3×N`

因为源码里更方便把每个原子存成一行：

$$
C = R^T
$$

所以同一个公式在源码里的行向量形式变成：

$$
C' = (C - \mathbf 1\mathbf r_0^T)(\mathrm{unit}\cdot A)^T
$$

这就是为什么你会看到 `axes` 出现在 `einsum('ix,kx->ki', ...)` 里，而不是直接写成你更熟悉的列向量左乘形式。

## 这和 `einsum` 怎么对应

源码里：

```python
c = numpy.einsum('ix,kx->ki', axes * unit, c - origin)
```

按下标写出来就是：

$$
c_{ki}^{\mathrm{new}} = \sum_x (A_{ix}\cdot \mathrm{unit})(c_{kx}-\mathrm{origin}_x)
$$

这里：

- `k`：第几个原子
- `x`：旧坐标分量
- `i`：新坐标分量

因此它正是：

```text
对每个原子 k
    把旧坐标在 x/y/z 三个方向上的分量
    投影到新坐标轴 i 上
```

如果你只看单个原子的列向量语言，它就是：

$$
\mathbf r_k' = (\mathrm{unit}\cdot A)(\mathbf r_k-\mathbf r_0)
$$

如果你看源码里的所有原子行向量形式，它就是：

```python
c = (c - origin) @ (axes * unit).T
```

这和 `einsum` 完全等价。

## MYNOTE：读 `format_atom` 时要抓住什么

看到 `format_atom` 这一段时，最重要的是抓住：

```text
origin: 定义新原点
axes:   定义新坐标轴在旧坐标系中的方向
unit:   把输入长度统一换成 Bohr
```

于是最后那一行并不是神秘操作，而只是：

```text
平移 + 换轴/旋转 + 单位变换
```