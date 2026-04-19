# symmetry_subgroup 概念笔记

## 一句话理解

`symmetry_subgroup` 里的 `subgroup` 是点群理论中的“子群”：它表示在分子的完整对称性中，只选择其中一部分对称操作来作为实际计算使用的点群。

换句话说：

> 分子本身可能有更高点群，但 PySCF 可以在计算时故意使用它的一个较低对称性的子群。

## topgroup 和 groupname

在 PySCF 中，和这个概念最相关的是两个属性：

| 属性 | 含义 |
|---|---|
| `mol.topgroup` | PySCF 从分子几何中检测到的最高点群 |
| `mol.groupname` | PySCF 实际用于计算的点群，可能是 `topgroup`，也可能是它的子群 |

例如，线性 CO 分子的真实对称性可能很高，比如 `Coov`。如果直接使用自动对称性：

```python
mol = gto.M(
    atom='C 0 0 0; O 0 0 1.5',
    symmetry=True,
)
```

PySCF 会检测最高点群，并选择一个实际可用的点群。

如果写成：

```python
mol = gto.M(
    atom='C 0 0 0; O 0 0 1.5',
    symmetry=True,
    symmetry_subgroup='C2v',
)
```

含义就是：

> 先让 PySCF 自动检测分子的最高对称性，但实际计算时使用 `C2v` 这个子群。

此时概念上可能出现：

```text
topgroup = Coov
groupname = C2v
```

## symmetry 和 symmetry_subgroup 的区别

| 写法 | 含义 |
|---|---|
| `symmetry=False` | 不使用点群对称性 |
| `symmetry=True` | 自动检测并使用对称性 |
| `symmetry=True, symmetry_subgroup='C2v'` | 自动检测最高点群，但实际使用指定子群 |
| `symmetry='C2v'` | 直接指定使用 `C2v` |

注意：如果 `symmetry` 已经被显式指定成字符串，例如 `symmetry='C2v'`，那么 `symmetry_subgroup` 就不再起作用。PySCF 的 `examples/gto/13-symmetry.py` 中也强调了这一点。

## 为什么需要子群

使用子群通常是为了程序实现和数值计算上的方便。

完整点群有时太高、太复杂，或者包含程序不方便直接处理的表示。例如线性分子可能对应 `Dooh` 或 `Coov`，但量子化学程序常常更方便使用有限的、阿贝尔的子群，例如：

```text
Dooh -> D2h
Coov -> C2v
```

原因包括：

- 阿贝尔群的不可约表示是一维的，轨道分块更简单。
- 程序可以更容易给轨道标记 irrep。
- 某些高对称群或无限群在具体积分和矩阵分块中实现更麻烦。
- 降低对称性有时能避免简并表示带来的处理复杂度。

## 在 mole.py 中的源码位置

在 `pyscf/gto/mole.py` 中，相关逻辑大致是：

```python
self.topgroup, orig, axes = symm.detect_symm(self._atom, self._basis)
```

这一步检测最高点群。

如果 `symmetry=True`，并且设置了 `symmetry_subgroup`，PySCF 会调用类似逻辑：

```python
groupname, axes = symm.as_subgroup(
    self.topgroup, axes, self.symmetry_subgroup
)
```

这一步把最高点群转换成指定子群。最后：

```python
self.groupname = groupname
```

也就是说，`topgroup` 是检测结果，`groupname` 是实际使用结果。


## 第一周需要掌握到什么程度

现在不需要深入学习所有点群操作，也不需要背子群表。第一周只要知道：

1. `symmetry_subgroup` 的 `sub` 指的是点群子群。
2. 子群意味着只使用完整对称性的一部分。
3. `topgroup` 是检测到的最高点群。
4. `groupname` 是实际计算使用的点群。
5. 如果显式写 `symmetry='C2v'`，`symmetry_subgroup` 通常不再起作用。