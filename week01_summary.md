# Week 01 Summary：从 Mole 输入到 libcint 数据结构

## 1. 本周主线

第一周的核心不是把 `pyscf/gto/mole.py` 全部读完，而是把下面这条数据流走通：

```text
用户输入 atom / basis / charge / spin / unit
    -> format_atom / format_basis
    -> self._atom / self._basis
    -> make_env
    -> self._atm / self._bas / self._env
    -> libcint 可以读取的底层数组
```

现在可以把 `Mole` 理解成一个“输入编译器”：它把用户熟悉的化学输入，逐步转换成积分库能够直接使用的数组结构。

---

## 2. 本周已经掌握的核心问题

### `mol.atom` 怎样变成 `self._atom`

`mol.atom` 是用户输入的原子和坐标。`format_atom()` 会做几件事：

```text
解析原子符号和坐标
处理输入单位
平移 origin
按 axes 做坐标变换
把内部坐标统一转成 Bohr
```

最后得到 `self._atom`，形如：

```python
[
    ['O', [0.0, 0.0, 0.0]],
    ['H', [0.0, 0.0, 1.8]],
]
```

这里的坐标已经是 PySCF 内部使用的 Bohr 单位。

相关笔记：

- `code_reading_notes/format_atom_coordinate_transform.md`
- `python_knowledge_notes/einsum.md`

---

### `mol.basis` 怎样变成 `self._basis`

`mol.basis` 可以是字符串、字典、按元素指定的 basis 等。构建过程中主要经过两步：

```text
_parse_default_basis
    -> 处理 default / 按元素分配这层外壳

format_basis
    -> 把 basis 数据解析成 PySCF 内部 shell 列表
```

`self._basis` 是一个字典：

```python
{
    'H': [[0, [alpha1, c1], [alpha2, c2], ...]],
    'O': [[0, ...], [0, ...], [1, ...], ...],
}
```

这里每个 value 是某个元素的一组 shell。后面的 `make_bas_env(basis_add, ...)` 中，`basis_add` 就是这个字典里某一个元素对应的 value。

相关笔记：

- `code_reading_notes/parse_default_basis.md`
- `concept_notes/basis_numbers_and_gto_formula.md`

---

### basis 数字怎样对应到 GTO 公式

PySCF 内部 basis shell 大致长这样：

```python
[l,
    [alpha_1, c_11, c_12, ...],
    [alpha_2, c_21, c_22, ...],
]
```

可以理解成：

```text
l        -> 角动量
alpha_p  -> primitive Gaussian exponent
c_pk     -> contraction coefficient
nprim    -> primitive 个数
nctr     -> contraction 个数，也就是 coefficient 矩阵列数
```

一个 shell 不是一个 AO。一个 shell 的 AO 数还要乘上角向函数个数。

例如 p shell 如果 `nctr=1`，它不是产生 1 个 AO，而是产生 3 个 p 型 AO。

相关笔记：

- `concept_notes/basis_numbers_and_gto_formula.md`
- `concept_notes/nao_and_ao_count.md`

---

### `_atm / _bas / _env` 各自是什么

这是第一周最重要的内部结构。

```text
_atm = 原子目录表
_bas = shell 目录表
_env = 浮点数据池
```

更具体地说：

```text
_atm 记录：原子核电荷、坐标在 _env 中的位置、核模型等
_bas 记录：shell 属于哪个原子、角动量、nprim、nctr、exponent/coeff 在 _env 中的位置
_env 记录：真正的浮点数，例如坐标、exponent、coefficient
```

`_atm` 和 `_bas` 像目录，`_env` 像仓库。

`libcint` 不直接吃 Python 的 list/dict/object，而是读取这些紧凑的数组。

相关笔记：

- `code_reading_notes/make_env_h2o_example.md`
- `python_knowledge_notes/array_stacking_and_concatenation.md`

---

### `make_env()` 做了什么

`make_env()` 可以拆成三步：

```text
1. 遍历 atoms
   -> make_atm_env
   -> 生成 _atm 行和坐标相关 _env 数据

2. 遍历 basis.items()
   -> make_bas_env
   -> 按元素生成 basis 模板和 exponent/coefficient 数据

3. 再次遍历 atoms
   -> 把 basis 模板复制给每个具体原子
   -> 修改 _bas 的 ATOM_OF 列
```

关键理解：

```text
make_bas_env 生成的是“元素模板”
最后那个 for 循环生成的是“原子实例”
```

所以两个 H 原子可以共用同一份 exponent / coefficient 数据，但它们在 `_bas` 中通过不同的 `ATOM_OF` 挂到不同核中心上。

---

### GTO 归一化怎么在 PySCF 中处理

`make_bas_env()` 中有两层归一化：

```python
cs = numpy.einsum('pi,p->pi', cs, gto_norm(angl, es))
if NORMALIZE_GTO:
    cs = _nomalize_contracted_ao(angl, es, cs)
```

可以理解成：

```text
原始 basis coefficient
    -> 乘上 primitive GTO normalization
    -> 再乘上每个 contracted function 的整体归一化因子
    -> 存进 _env 给 libcint 直接使用
```

`bas_ctr_coeff()` 会把 primitive normalization 除掉，但仍保留 contracted normalization。

相关笔记：

- `concept_notes/gto_normalization.md`

---

### `charge / spin / nelectron / nelec / multiplicity` 的关系

PySCF 的 `spin` 不是多重度，而是：

$$
\mathrm{spin} = N_\alpha - N_\beta = 2S
$$

多重度是：

$$
\mathrm{multiplicity} = 2S + 1 = \mathrm{spin} + 1
$$

总电子数是：

$$
N_e = \sum_A Z_A^{\mathrm{eff}} - q
$$

其中 $q$ 是分子净电荷，对应 `mol.charge`。

`build()` 中的：

```python
self.nelec
```

虽然没有接收返回值，但它会触发 `nelec` property 的 getter，用来检查 `nelectron` 和 `spin` 是否自洽。

相关笔记：

- `concept_notes/electron_count_spin_multiplicity.md`

---

### `nao` 怎样从 shell 算出来

`nao` 在 `mol.nao_nr()` 这类接口中主要表示：

```text
number of atomic orbitals
```

非相对论 spherical 情况下：

$$
N_{\mathrm{AO}}
=
\sum_{\mathrm{shells}} (2l+1)n_{\mathrm{ctr}}
$$

Cartesian 情况下：

$$
N_{\mathrm{AO}}
=
\sum_{\mathrm{shells}} \frac{(l+1)(l+2)}{2}n_{\mathrm{ctr}}
$$

所以：

```text
nbas = shell 数
nprim = primitive 个数
nctr = contracted radial function 个数
nao = AO basis function 个数
```

这几个量不能混用。

相关笔记：

- `concept_notes/nao_and_ao_count.md`

---

## 3. 本周涉及的 Python / NumPy 工具

### `exec`

用于执行字符串形式的 Python 代码。PySCF 的 `__config__.py` 会用它读取用户配置文件 `pyscf_conf.py`，从而覆盖默认配置值。

### `tempfile`

Python 标准库模块，用来创建临时文件或临时目录。

### `property`

`@property` 可以让方法像属性一样被访问。例如：

```python
self.nelec
```

实际上会执行：

```python
nelec(self)
```

这解释了为什么 `build()` 中单独访问 `self.nelec` 能触发一致性检查。

### `einsum`

`numpy.einsum()` 用下标记号表达求和和矩阵运算。本周主要遇到：

```python
numpy.einsum('ix,kx->ki', axes * unit, c - origin)
numpy.einsum('pi,p->pi', cs, gto_norm(angl, es))
numpy.einsum('pi,pq,qi->i', cs, ee, cs)
```

### `vstack / hstack`

在 `make_env()` 中：

```text
vstack 用来把 _atm / _bas 的行记录堆成二维表
hstack 用来把 _env 的许多浮点片段接成一维数据池
```

相关笔记：

- `python_knowledge_notes/property_setter.md`
- `python_knowledge_notes/einsum.md`
- `python_knowledge_notes/array_stacking_and_concatenation.md`

---

## 4. 本周提前接触但不必继续深挖的内容

这些内容已经有基本笔记，但可以先不作为第二周主线：

```text
symmetry_subgroup
irreducible representation
symmetry-adapted AO basis
NAO 的另一个含义：Natural Atomic Orbital
ECP / pseudo 的细节
spinor / nao_2c
```

相关笔记：

- `concept_notes/symmetry_subgroup.md`
- `concept_notes/symmetry_adapted_ao_basis.md`

---

## 5. 第一周完成标准

现在应该能够回答这些问题：

1. `mol.atom` 怎样变成 `self._atom`？
2. `mol.basis` 怎样变成 `self._basis`？
3. `_atm / _bas / _env` 各自存什么？
4. `make_env()` 为什么要先做 basis 模板，再分配给具体原子？
5. basis 中的 exponent 和 coefficient 怎样对应到 GTO 公式？
6. primitive normalization 和 contracted normalization 分别是什么？
7. `charge / spin / nelec / multiplicity` 怎样互相约束？
8. 一个 shell 怎样贡献 AO 数？
9. H2O 中 O 用 6-31G、H 用 STO-3G 时，为什么 `nao = 11`？
10. 为什么 PySCF 最后要把数据打包成 `_atm / _bas / _env` 给 `libcint`？

如果这些问题可以用自己的话讲出来，第一周就已经完成。

---

## 6. 第二周入口

第二周最自然的主题是：

```text
AO integrals：PySCF 怎样从 _atm/_bas/_env 生成积分矩阵
```

建议从这些问题开始：

```text
mol.intor() 是什么？
int1e_ovlp / int1e_kin / int1e_nuc 分别对应什么积分？
积分矩阵为什么是 (nao, nao)？
shls_slice 如何按 shell 截取积分？
PySCF 如何把 intor 名字传给 libcint？
```

第二周不需要马上读 SCF 主循环。先把 AO 积分这一层接上，会让后面读 RHF / RKS 更稳。

---

## 7. 一句话收尾

第一周真正建立起来的能力是：

```text
看到一个 PySCF Mole 输入，能够解释它怎样一步步变成 libcint 使用的底层数组，并能从 shell 结构推断 AO 数和矩阵维度。
```