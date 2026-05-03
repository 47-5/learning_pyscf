# NumPy / PyTorch 中的数组拼接与堆叠：vstack、hstack、concatenate、stack

## 1. 先记住一个核心问题

看到 `vstack`、`hstack`、`concatenate`、`stack` 这类函数时，先不要急着背名字，先问自己一句：

```text
我是沿着已有的轴继续接，还是新建一根轴来装这些数组？
```

这句话基本可以把最常见的混乱分成两类：

```text
concatenate / hstack / vstack / cat  -> 沿已有轴拼接
stack                                -> 新建一根轴
```

对二维数组来说，习惯上：

```text
axis=0  行方向，往下接，行数增加
axis=1  列方向，往右接，列数增加
```

例如一个矩阵形状是 `(2, 3)`：

```text
axis=0:  有 2 行
axis=1:  有 3 列
```

如果沿 `axis=0` 拼接，结果会变成更多行；如果沿 `axis=1` 拼接，结果会变成更多列。

---

## 2. `np.concatenate`：沿已有轴拼接

`np.concatenate` 是最基础、最值得先理解的函数。

```python
import numpy as np

a = np.array([[1, 2],
              [3, 4]])

b = np.array([[5, 6],
              [7, 8]])
```

### 沿 `axis=0` 拼接

```python
np.concatenate([a, b], axis=0)
```

结果是：

```python
array([[1, 2],
       [3, 4],
       [5, 6],
       [7, 8]])
```

形状变化：

```text
(2, 2) + (2, 2) -> (4, 2)
```

行数增加，列数不变。

### 沿 `axis=1` 拼接

```python
np.concatenate([a, b], axis=1)
```

结果是：

```python
array([[1, 2, 5, 6],
       [3, 4, 7, 8]])
```

形状变化：

```text
(2, 2) + (2, 2) -> (2, 4)
```

列数增加，行数不变。

### 形状要求

沿某个轴拼接时，除了那个轴可以不同，其他轴必须相同。

比如沿 `axis=0` 拼接时：

```text
(2, 3) 和 (5, 3) 可以拼 -> (7, 3)
(2, 3) 和 (5, 4) 不可以拼，因为列数不同
```

沿 `axis=1` 拼接时：

```text
(2, 3) 和 (2, 5) 可以拼 -> (2, 8)
(2, 3) 和 (4, 5) 不可以拼，因为行数不同
```

---

## 3. `np.vstack`：按行往下堆

`vstack` 可以理解成 vertical stack，也就是竖直方向堆叠。

对二维数组，它基本等价于：

```python
np.concatenate([a, b], axis=0)
```

例如：

```python
a = np.array([[1, 2],
              [3, 4]])

b = np.array([[5, 6]])

np.vstack([a, b])
```

结果：

```python
array([[1, 2],
       [3, 4],
       [5, 6]])
```

形状变化：

```text
(2, 2) + (1, 2) -> (3, 2)
```

### 和 PySCF 中 `_atm` / `_bas` 的关系

在 `make_env()` 的最后，PySCF 有类似这样的代码：

```python
_atm = numpy.asarray(numpy.vstack(_atm), numpy.int32).reshape(-1, ATM_SLOTS)
_bas = numpy.asarray(numpy.vstack(_bas), numpy.int32).reshape(-1, BAS_SLOTS)
```

这很好理解：

```text
_atm 是很多个“原子记录行”堆起来的表
_bas 是很多个“shell 记录行”堆起来的表
```

所以这里用 `vstack` 很自然，因为它想把许多行记录往下堆成二维表。

假设：

```python
atm0 = np.array([1, 20, 1, 23, 0, 0])
atm1 = np.array([1, 24, 1, 27, 0, 0])
```

那么：

```python
np.vstack([atm0, atm1])
```

得到：

```python
array([[ 1, 20,  1, 23,  0,  0],
       [ 1, 24,  1, 27,  0,  0]])
```

这就是 `_atm` 表的雏形。

### `vstack` 对一维数组的特殊处理

这是很容易忘的一点。

```python
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

np.vstack([a, b])
```

结果不是一维数组，而是二维数组：

```python
array([[1, 2, 3],
       [4, 5, 6]])
```

形状变化：

```text
(3,) + (3,) -> (2, 3)
```

也就是说，`vstack` 会先把一维数组当成一行：

```text
(3,) -> (1, 3)
```

然后再往下堆。

---

## 4. `np.hstack`：按列往右接

`hstack` 可以理解成 horizontal stack，也就是水平方向拼接。

对二维数组，它基本等价于：

```python
np.concatenate([a, b], axis=1)
```

例如：

```python
a = np.array([[1, 2],
              [3, 4]])

b = np.array([[5],
              [6]])

np.hstack([a, b])
```

结果：

```python
array([[1, 2, 5],
       [3, 4, 6]])
```

形状变化：

```text
(2, 2) + (2, 1) -> (2, 3)
```

### 和 PySCF 中 `_env` 的关系

在 `make_env()` 最后，PySCF 有类似这样的代码：

```python
_env = numpy.asarray(numpy.hstack(_env), dtype=numpy.float64)
```

这里 `_env` 本来是一个列表，里面装了很多段浮点数据，例如：

```text
pre_env
原子 0 的坐标和 zeta
原子 1 的坐标和 zeta
某个元素的 exponents
某个元素的 coefficients
...
```

这些东西最后要变成一个一维大数组，也就是 libcint 使用的浮点数据池。

所以这里用 `hstack` 的直觉是：

```text
把很多段一维数据，首尾相接，接成一条长的一维 _env
```

例如：

```python
env0 = np.array([0.0, 0.0, 0.0, 0.0])
env1 = np.array([0.0, 0.0, 1.4, 0.0])
exp  = np.array([3.425, 0.624, 0.169])

np.hstack([env0, env1, exp])
```

结果：

```python
array([0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 1.4, 0.0,
       3.425, 0.624, 0.169])
```

### `hstack` 对一维数组的特殊处理

对一维数组，`hstack` 会直接拼成更长的一维数组：

```python
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

np.hstack([a, b])
```

结果：

```python
array([1, 2, 3, 4, 5, 6])
```

形状变化：

```text
(3,) + (3,) -> (6,)
```

这和 `vstack` 对一维数组的行为不同：

```python
np.vstack([a, b]).shape  # (2, 3)
np.hstack([a, b]).shape  # (6,)
```

这个差别很重要。

---

## 5. `np.stack`：新建一根轴

`stack` 和 `concatenate` 最大的区别是：

```text
concatenate 沿已有轴拼接
stack 新建一根轴
```

还是用两个一维数组：

```python
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])
```

### `axis=0`

```python
np.stack([a, b], axis=0)
```

结果：

```python
array([[1, 2, 3],
       [4, 5, 6]])
```

形状变化：

```text
(3,) 和 (3,) -> (2, 3)
```

这里新建的轴放在最前面，表示“第几个数组”。

### `axis=1`

```python
np.stack([a, b], axis=1)
```

结果：

```python
array([[1, 4],
       [2, 5],
       [3, 6]])
```

形状变化：

```text
(3,) 和 (3,) -> (3, 2)
```

这里新建的轴放在后面，表示“每个位置上收集两个数组的值”。

### 和 `vstack` 的相似与不同

对于两个一维数组：

```python
np.vstack([a, b])
np.stack([a, b], axis=0)
```

结果看起来一样，都是：

```python
array([[1, 2, 3],
       [4, 5, 6]])
```

但概念不同：

```text
vstack 是先把一维数组看成行，再沿 axis=0 拼接
stack 是直接新建一根 axis=0
```

对高维数组时，这个差别会更明显。

---

## 6. `column_stack` 和 `row_stack`

这两个函数不一定每天用，但看代码时偶尔会遇到。

### `np.column_stack`

把一维数组当成列来拼：

```python
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

np.column_stack([a, b])
```

结果：

```python
array([[1, 4],
       [2, 5],
       [3, 6]])
```

形状变化：

```text
(3,) + (3,) -> (3, 2)
```

它很适合“我有几列数据，想组成一张表”的场景。

### `np.row_stack`

`row_stack` 基本就是 `vstack` 的别名思路：把东西按行堆起来。

```python
np.row_stack([a, b])
```

结果：

```python
array([[1, 2, 3],
       [4, 5, 6]])
```

---

## 7. `dstack`：沿第三维堆叠

`dstack` 是 depth stack，可以先粗略理解成沿第三个轴拼。

```python
a = np.array([[1, 2],
              [3, 4]])

b = np.array([[5, 6],
              [7, 8]])

np.dstack([a, b])
```

结果形状是：

```text
(2, 2, 2)
```

你可以把它想成：

```text
每个位置 (i, j) 上，存了来自 a 和 b 的两个值
```

这在图像数据里比较常见，比如把不同 channel 拼到一起。

对读 PySCF 的 `mole.py` 来说，`dstack` 暂时不是主线，知道它是“第三维方向的 stack”就够了。

---

## 8. NumPy 小抄

| 函数 | 核心含义 | 常见等价理解 | 是否新建轴 |
|---|---|---|---|
| `np.concatenate(arrs, axis=0)` | 沿已有 `axis=0` 拼 | 往下接 | 否 |
| `np.concatenate(arrs, axis=1)` | 沿已有 `axis=1` 拼 | 往右接 | 否 |
| `np.vstack(arrs)` | 竖直堆叠 | 2D 时近似 `axis=0` 拼 | 否 |
| `np.hstack(arrs)` | 水平拼接 | 2D 时近似 `axis=1` 拼 | 否 |
| `np.stack(arrs, axis=0)` | 新建第 0 轴 | 多包一层 | 是 |
| `np.stack(arrs, axis=1)` | 新建第 1 轴 | 插入一根新轴 | 是 |
| `np.column_stack(arrs)` | 一维数组按列组成表 | 每个输入是一列 | 会把 1D 转成列 |
| `np.dstack(arrs)` | 沿第三维堆叠 | 深度/channel 方向 | 视输入维度而定 |

最值得记的不是函数名，而是这两句：

```text
concatenate / vstack / hstack：数组维度数量通常不变，只是某个轴变长
stack：数组维度数量会增加 1
```

---

## 9. PyTorch 中怎么对应

PyTorch 的核心思路和 NumPy 非常像，只是对象从 `ndarray` 换成了 `Tensor`。

### `torch.cat`：对应 NumPy 的 `concatenate`

```python
import torch

x = torch.tensor([[1, 2],
                  [3, 4]])

y = torch.tensor([[5, 6],
                  [7, 8]])
```

沿 `dim=0` 拼：

```python
torch.cat([x, y], dim=0)
```

结果：

```python
tensor([[1, 2],
        [3, 4],
        [5, 6],
        [7, 8]])
```

沿 `dim=1` 拼：

```python
torch.cat([x, y], dim=1)
```

结果：

```python
tensor([[1, 2, 5, 6],
        [3, 4, 7, 8]])
```

所以：

```text
np.concatenate(..., axis=0)  <->  torch.cat(..., dim=0)
np.concatenate(..., axis=1)  <->  torch.cat(..., dim=1)
```

在 PyTorch 里，`torch.concat` 和 `torch.concatenate` 也可以看成 `torch.cat` 的别名。

### `torch.stack`：新建一根维度

```python
a = torch.tensor([1, 2, 3])
b = torch.tensor([4, 5, 6])
```

```python
torch.stack([a, b], dim=0)
```

结果：

```python
tensor([[1, 2, 3],
        [4, 5, 6]])
```

形状：

```text
(3,) 和 (3,) -> (2, 3)
```

```python
torch.stack([a, b], dim=1)
```

结果：

```python
tensor([[1, 4],
        [2, 5],
        [3, 6]])
```

形状：

```text
(3,) 和 (3,) -> (3, 2)
```

这和 NumPy 的 `np.stack` 思路一样。

### `torch.vstack` 和 `torch.hstack`

PyTorch 也有 `torch.vstack` 和 `torch.hstack`。

```python
a = torch.tensor([1, 2, 3])
b = torch.tensor([4, 5, 6])
```

```python
torch.vstack([a, b])
```

结果：

```python
tensor([[1, 2, 3],
        [4, 5, 6]])
```

```python
torch.hstack([a, b])
```

结果：

```python
tensor([1, 2, 3, 4, 5, 6])
```

对二维张量：

```python
x = torch.tensor([[1],
                  [2],
                  [3]])

y = torch.tensor([[4],
                  [5],
                  [6]])
```

```python
torch.hstack([x, y])
```

结果：

```python
tensor([[1, 4],
        [2, 5],
        [3, 6]])
```

这个行为和 NumPy 很接近：

```text
1D: hstack 直接接成长向量
2D: hstack 按列往右接
```

---

## 10. PyTorch 小抄

| NumPy | PyTorch | 说明 |
|---|---|---|
| `np.concatenate(arrs, axis=0)` | `torch.cat(tensors, dim=0)` | 沿已有第 0 维拼接 |
| `np.concatenate(arrs, axis=1)` | `torch.cat(tensors, dim=1)` | 沿已有第 1 维拼接 |
| `np.stack(arrs, axis=0)` | `torch.stack(tensors, dim=0)` | 新建第 0 维 |
| `np.stack(arrs, axis=1)` | `torch.stack(tensors, dim=1)` | 新建第 1 维 |
| `np.vstack(arrs)` | `torch.vstack(tensors)` | 竖直堆叠 |
| `np.hstack(arrs)` | `torch.hstack(tensors)` | 水平拼接 |

最重要的对应关系是：

```text
NumPy 的 axis 约等于 PyTorch 的 dim
```

---

## 11. `cat` 和 `stack` 的区别：最容易考自己的例子

假设：

```python
x.shape == (2, 3)
y.shape == (2, 3)
```

### 沿已有轴拼接

```python
np.concatenate([x, y], axis=0).shape
# (4, 3)

np.concatenate([x, y], axis=1).shape
# (2, 6)
```

PyTorch 对应：

```python
torch.cat([x, y], dim=0).shape
# torch.Size([4, 3])

torch.cat([x, y], dim=1).shape
# torch.Size([2, 6])
```

### 新建轴

```python
np.stack([x, y], axis=0).shape
# (2, 2, 3)

np.stack([x, y], axis=1).shape
# (2, 2, 3)

np.stack([x, y], axis=2).shape
# (2, 3, 2)
```

注意，前两个例子形状都写成 `(2, 2, 3)`，但含义不同：

```text
axis=0: 新轴放在最前面，表示第几个数组
axis=1: 新轴插在原来的第 0 轴和第 1 轴之间
```

对简单形状来说它们可能碰巧一样，但轴的语义不一样。真正写代码时，后续索引方式会不同。

---

## 12. 回到 PySCF：为什么 `_atm/_bas` 用 `vstack`，`_env` 用 `hstack`

在 `make_env()` 结尾，PySCF 做的是两种不同的数据整理。

### `_atm` 和 `_bas` 是二维表

```text
_atm: 每一行是一个原子记录
_bas: 每一行是一个 shell 记录
```

所以它们要往下堆成表：

```python
numpy.vstack(_atm)
numpy.vstack(_bas)
```

直觉：

```text
一行
一行
一行
往下堆
```

### `_env` 是一维数据池

```text
_env: 一段坐标 + 一段 zeta + 一段 exponent + 一段 coefficient + ...
```

它不是表，而是一条长长的浮点数组。所以要把很多段首尾相接：

```python
numpy.hstack(_env)
```

直觉：

```text
一段 + 一段 + 一段 + 一段 -> 一条长数组
```

这就是为什么同一个函数结尾会同时出现：

```python
vstack: 整理二维记录表
hstack: 整理一维数据池
```

---

## 13. 读代码时的判断流程

以后看到类似代码，可以按这个顺序判断：

1. 先看每个元素的 shape
2. 再问结果应该是“表”还是“一条长数组”
3. 如果是表，通常考虑 `vstack` 或 `concatenate(axis=0)`
4. 如果是横向扩列，考虑 `hstack` 或 `concatenate(axis=1)`
5. 如果需要多出一个“第几个样本/第几个数组”的维度，考虑 `stack`

对你当前读 PySCF 的主线，最重要的是：

```text
vstack 让很多行记录变成二维表
hstack 让很多段数据变成一维数据池
stack 会新建轴，当前 make_env 这段并不需要它
```

---

## 14. 参考文档

- NumPy `concatenate`: https://numpy.org/doc/stable/reference/generated/numpy.concatenate.html
- NumPy `vstack`: https://numpy.org/doc/stable/reference/generated/numpy.vstack.html
- NumPy `hstack`: https://numpy.org/doc/stable/reference/generated/numpy.hstack.html
- NumPy `stack`: https://numpy.org/doc/stable/reference/generated/numpy.stack.html
- PyTorch `torch.cat`: https://docs.pytorch.org/docs/stable/generated/torch.cat.html
- PyTorch `torch.stack`: https://docs.pytorch.org/docs/stable/generated/torch.stack.html
- PyTorch `torch.vstack`: https://docs.pytorch.org/docs/stable/generated/torch.vstack.html
- PyTorch `torch.hstack`: https://docs.pytorch.org/docs/stable/generated/torch.hstack.html