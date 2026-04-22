# `_parse_default_basis` 与 basis / ECP / pseudo 输入管线

## 一句话理解

`_parse_default_basis` 的名字有一点误导。它并不真正解析 Gaussian basis、ECP 或 pseudo 的具体数据格式。

它只做一件事：

> 把“所有原子共用一个默认设置”的输入，统一展开成 `{atom: spec}` 这种按元素分发的字典形式。

真正把数据变成内部格式的是后面的三个函数：

```python
format_basis(...)
format_ecp(...)
format_pseudo(...)
```

所以：

```text
basis / ecp / pseudo
    共享外层输入分发模式
    但内部数据格式完全不同
```

## 在 `Mole.build()` 里的三条路径

在 `pyscf/gto/mole.py` 的 `Mole.build()` 中，可以看到三条相似但不相同的路径。

### basis 路径

```python
if self.basis:
    _basis = _parse_default_basis(self.basis, uniq_atoms)
    self._basis = self.format_basis(_basis)
```

含义：

```text
用户输入 self.basis
    -> 按元素展开为 {atom: basis_spec}
    -> format_basis
    -> self._basis
```

`self._basis` 是 Python 层内部使用的 AO basis 格式，后续会被 `make_env()` 转成 libcint 需要的 `_bas` 和 `_env`。

### ECP 路径

```python
if self.ecp:
    atoms_wo_ghost = [a for a in uniq_atoms if not is_ghost_atom(a)]
    _ecp = _parse_default_basis(self.ecp, atoms_wo_ghost)
    self._ecp = self.format_ecp(_ecp)
```

含义：

```text
用户输入 self.ecp
    -> 按元素展开为 {atom: ecp_spec}
    -> format_ecp
    -> self._ecp
```

`self._ecp` 描述的是 effective core potential，也就是用有效势替代芯电子的那部分信息。

### pseudo 路径

```python
if self.pseudo:
    atoms_wo_ghost = [a for a in uniq_atoms if not is_ghost_atom(a)]
    _pseudo = _parse_default_basis(self.pseudo, atoms_wo_ghost)
    self._pseudo = _pseudo = self.format_pseudo(_pseudo)
```

含义：

```text
用户输入 self.pseudo
    -> 按元素展开为 {atom: pseudo_spec}
    -> format_pseudo
    -> self._pseudo
```

`self._pseudo` 描述的是 pseudopotential，尤其常见于 PBC / GTH / CP2K 风格的周期体系计算。

## `_parse_default_basis` 做了什么

源码大致是：

```python
def _parse_default_basis(basis, uniq_atoms):
    if isinstance(basis, (str, tuple, list)):
        _basis = {a: basis for a in uniq_atoms}
    elif 'default' in basis:
        default_basis = basis['default']
        _basis = {a: default_basis for a in uniq_atoms}
        _basis.update(basis)
        del _basis['default']
    else:
        _basis = basis
    return _basis
```

这里的参数名虽然叫 `basis`，但它也会被传入 `self.ecp` 和 `self.pseudo`。

它支持三种外层写法。

### 写法 1：所有元素共用一个设置

```python
basis = 'sto-3g'
uniq_atoms = {'O', 'H'}
```

会变成：

```python
{
    'O': 'sto-3g',
    'H': 'sto-3g',
}
```

类似地：

```python
ecp = 'lanl2dz'
pseudo = 'gth-pbe'
```

也可以先被展开成按元素分发的字典。

### 写法 2：`default` 加局部覆盖

```python
basis = {
    'default': 'def2-svp',
    'O': 'cc-pvdz',
}
```

如果体系中有 O 和 H，会变成：

```python
{
    'O': 'cc-pvdz',
    'H': 'def2-svp',
}
```

也就是先给所有元素填默认值，再用显式元素设置覆盖。

### 写法 3：已经是按元素指定的字典

```python
basis = {
    'O': 'cc-pvdz',
    'H': 'sto-3g',
}
```

则直接返回原字典。

## 三者哪里开始不同

`_parse_default_basis` 之后，三者立刻分流。

| 输入属性 | 外层展开后 | 后续函数 | 内部含义 |
|---|---|---|---|
| `basis` | `{atom: basis_spec}` | `format_basis` | AO Gaussian basis shells |
| `ecp` | `{atom: ecp_spec}` | `format_ecp` | effective core potential |
| `pseudo` | `{atom: pseudo_spec}` | `format_pseudo` | pseudopotential, often GTH / PBC |

所以它们只是共享 `{atom: spec}` 这个外壳，不共享内部格式。

## `format_basis` 的结果

`format_basis()` 会把 AO basis 输入转换成类似：

```python
{
    'H': [
        [0,
         [3.42525091, 0.15432897],
         [0.62391373, 0.53532814],
         [0.16885540, 0.44463454]]
    ]
}
```

可以理解为：

```text
atom -> shell list
shell -> [angular_momentum, primitive1, primitive2, ...]
primitive -> [exponent, contraction_coefficients...]
```

这里的重点是：

> `format_basis` 处理的是 AO basis，也就是后续积分中真正定义 AO 函数空间的内容。

## `format_ecp` 的结果

`format_ecp()` 处理的是 effective core potential。它的内部结构大致是：

```text
atom -> (
    nelec,
    angular-momentum-dependent potential terms
)
```

其中 `nelec` 是被 ECP 替代掉的 core electron 数量。

这和 AO basis 完全不是同一种对象。AO basis 定义的是展开波函数用的函数空间；ECP 定义的是芯电子被替代后，价电子感受到的有效势。

## `format_pseudo` 的结果

`format_pseudo()` 处理 pseudopotential，尤其是周期体系中常见的 GTH / CP2K 风格 pseudo。

它的数据结构包含：

```text
atom -> local part + nonlocal projectors + effective valence charges
```

细节比第一周需要掌握的内容复杂。现在只要知道：

> `pseudo` 和 `ecp` 都是用势函数替代部分核心电子或核-电子相互作用的技术，但 PySCF 中它们走不同的数据格式和不同的后续处理路径。

## 该读哪里

第一周建议阅读顺序：

1. `pyscf/gto/mole.py` 中 `Mole.build()` 的三段调用。
2. `pyscf/gto/mole.py` 中 `_parse_default_basis`。
3. `pyscf/gto/mole.py` 中 `format_basis`。
4. 粗读 `format_ecp` 和 `format_pseudo` 的 docstring。
5. 暂时跳过 ECP/pseudo parser 的细节。

如果要追 loader：

| 目标 | 入口 |
|---|---|
| basis loader | `pyscf/gto/basis/__init__.py::load` |
| ECP loader | `pyscf/gto/basis/__init__.py::load_ecp` |
| pseudo loader | `pyscf/gto/basis/__init__.py::load_pseudo` |

注意：`pseudo` 的入口在 `pyscf/pbc/gto/pseudo/__init__.py` 中转接到 `pyscf.gto.basis.load_pseudo`。所以不要被 `basis` 这个包名迷惑。

## 可以写进源码注释的 MYNOTE

```python
# MYNOTE: _parse_default_basis does not parse basis/ECP/pseudo data itself.
# MYNOTE: It only expands a global/default specification into {atom: spec}.
# MYNOTE: basis, ecp, and pseudo share this outer per-atom dispatch pattern.
# MYNOTE: Their actual internal formats diverge in format_basis/format_ecp/format_pseudo.
```

## 第一周需要掌握到什么程度

现在只需要掌握：

1. `_parse_default_basis` 是外层分发工具，不是真正的 basis parser。
2. `basis`、`ecp`、`pseudo` 都可以用统一的按元素指定方式。
3. 三者经过 `_parse_default_basis` 后都变成 `{atom: spec}`。
4. 三者随后进入不同的 `format_*` 函数。
5. 第一周重点读 `_basis`，因为它直接决定 AO 空间和后续积分。
6. ECP/pseudo 先知道它们和 basis 共享外层模式，但内部格式和物理含义不同即可。