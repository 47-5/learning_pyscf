# `nao`：PySCF 中 AO 数、shell 数和 basis function 数的关系

## 1. 这里的 `nao` 是什么意思

在 PySCF 的 `mole.py` 里，`nao` 通常可以读成：

```text
number of atomic orbitals
```

也就是 AO basis functions 的数量。

这点要和另一个常见术语区分开：在量子化学文献里，`NAO` 有时也指 **Natural Atomic Orbital**。但在 `mol.nao_nr()`、`mol.nao`、`nao_cart()` 这些 PySCF 接口里，`nao` 主要是“AO 个数”的意思。

所以在读 `mole.py` 这一段时，可以先把它理解成：

```text
nao = 当前 Mole 对象中 AO basis function 的总数
```

---

## 2. 从 shell 到 AO 数的核心公式

PySCF 的 `_bas` 是按 shell 存的。每个 shell 至少有这些信息：

```text
ANG_OF   -> 角动量 l
NPRIM_OF -> primitive Gaussian 个数 nprim
NCTR_OF  -> contraction 个数 nctr
```

一个 shell 对 AO 数的贡献不是 `nprim`，而是：

```text
角向函数个数 × contraction 个数
```

### spherical GTO

对于 spherical harmonic AO，角动量为 $l$ 的 shell 有：

$$
d_l^{\mathrm{sph}} = 2l + 1
$$

个角向函数。

所以一个 shell 的 AO 数是：

$$
N_{\mathrm{AO,shell}}^{\mathrm{sph}}
= (2l+1)n_{\mathrm{ctr}}
$$

例如：

```text
s shell: l=0 -> 1 个角向函数
p shell: l=1 -> 3 个角向函数
d shell: l=2 -> 5 个角向函数
f shell: l=3 -> 7 个角向函数
```

### Cartesian GTO

对于 Cartesian GTO，角动量为 $l$ 的 shell 有：

$$
d_l^{\mathrm{cart}}
= \frac{(l+1)(l+2)}{2}
$$

个角向函数。

所以一个 shell 的 AO 数是：

$$
N_{\mathrm{AO,shell}}^{\mathrm{cart}}
= \frac{(l+1)(l+2)}{2}n_{\mathrm{ctr}}
$$

例如：

```text
s shell: l=0 -> 1 个 Cartesian 函数
p shell: l=1 -> 3 个 Cartesian 函数
d shell: l=2 -> 6 个 Cartesian 函数
f shell: l=3 -> 10 个 Cartesian 函数
```

因此 spherical 和 Cartesian 在 s、p shell 上数量相同，但从 d shell 开始不同：

```text
s: 1 vs 1
p: 3 vs 3
d: 5 vs 6
f: 7 vs 10
```

---

## 3. `nao_nr()`：非相对论 contracted AO 总数

源码里 `nao_nr()` 是：

```python
def nao_nr(mol, cart=None):
    '''Total number of contracted GTOs for the given Mole object'''
    if cart is None:
        cart = mol.cart
    if cart:
        return nao_cart(mol)
    else:
        return int(((mol._bas[:,ANG_OF]*2+1) * mol._bas[:,NCTR_OF]).sum())
```

这里的 `nr` 可以理解成 **non-relativistic**。

如果不是 Cartesian，它计算的是 spherical contracted GTO 个数：

$$
N_{\mathrm{AO}}^{\mathrm{sph}}
=
\sum_{\mathrm{shells}} (2l+1)n_{\mathrm{ctr}}
$$

也就是说，`nao_nr()` 默认数的是：

```text
非相对论、contracted、通常为 spherical 的 AO 数
```

如果 `cart=True`，它会转去调用 `nao_cart()`。

---

## 4. `nao_cart()`：Cartesian contracted AO 总数

源码里 `nao_cart()` 是：

```python
def nao_cart(mol):
    '''Total number of contracted cartesian GTOs for the given Mole object'''
    l = mol._bas[:,ANG_OF]
    return int(((l+1)*(l+2)//2 * mol._bas[:,NCTR_OF]).sum())
```

对应公式：

$$
N_{\mathrm{AO}}^{\mathrm{cart}}
=
\sum_{\mathrm{shells}}
\frac{(l+1)(l+2)}{2}n_{\mathrm{ctr}}
$$

它和 `nao_nr(cart=False)` 的区别主要体现在 d、f、g 等高角动量 shell。

---

## 5. `mol.nao` property：默认就是 `mol.nao_nr()`

`Mole` 类中还有一个 property：

```python
@property
def nao(self):
    if self._nao is None:
        return self.nao_nr()
    else:
        return self._nao

@nao.setter
def nao(self, x):
    self._nao = x
```

所以通常情况下：

```python
mol.nao == mol.nao_nr()
```

除非用户或某些特殊流程显式设置过：

```python
mol.nao = something
```

对于正常读 `mole.py`，可以先把 `mol.nao` 理解成：

```text
当前 Mole 对象的默认 AO 个数
```

---

## 6. `npgto_nr()`：primitive GTO 数，不是 AO 数

源码里还有一个名字很像的函数：

```python
def npgto_nr(mol, cart=None):
    '''Total number of primitive spherical GTOs for the given Mole object'''
    if cart is None:
        cart = mol.cart
    l = mol._bas[:,ANG_OF]
    if cart:
        return int(((l+1)*(l+2)//2 * mol._bas[:,NPRIM_OF]).sum())
    else:
        return int(((l*2+1) * mol._bas[:,NPRIM_OF]).sum())
```

它数的是 primitive GTO 数，因此用的是 `NPRIM_OF`，不是 `NCTR_OF`。

对比一下：

```text
nao_nr   -> 用 NCTR_OF，数 contracted AO
npgto_nr -> 用 NPRIM_OF，数 primitive GTO
```

这和我们前面 basis 笔记里的区分一致：

```text
primitive 个数决定原始 Gaussian 数量
contraction 个数决定 contracted radial function 数量
AO 数还要再乘角向函数个数
```

---

## 7. 一个 shell 的贡献如何手算

假设某个 p shell：

```text
l = 1
nprim = 3
nctr = 1
```

如果是 spherical 或 Cartesian，p shell 的角向函数个数都是 3。

所以：

$$
N_{\mathrm{AO,shell}} = 3 \times 1 = 3
$$

这三个 AO 可以理解成：

```text
p_x, p_y, p_z
```

或者对应的三个 real spherical p functions。

注意，`nprim=3` 不表示这个 shell 产生 3 个 AO。它只是说这个 contracted radial function 是由 3 个 primitive Gaussian 收缩出来的。

再看一个 d shell：

```text
l = 2
nprim = 6
nctr = 1
```

spherical 下：

$$
N_{\mathrm{AO,shell}}^{\mathrm{sph}} = (2\times2+1)\times1=5
$$

Cartesian 下：

$$
N_{\mathrm{AO,shell}}^{\mathrm{cart}} = \frac{(2+1)(2+2)}{2}\times1=6
$$

这就是常说的 5D / 6D 差别。

---

## 8. H2O 例子：O 用 6-31G，H 用 STO-3G

我们前面用过一个混合 basis 的 H2O 例子：

```python
basis = {
    'O': '6-31g',
    'H': 'sto-3g',
}
```

简化地看，O 的 6-31G 可以数成：

```text
O: [3s, 2p]
```

也就是：

```text
3 个 s 型 contracted radial functions
2 组 p 型 contracted radial functions
```

对于 spherical AO：

```text
O 的 s 部分: 3 × 1 = 3
O 的 p 部分: 2 × 3 = 6
O 总 AO 数:  3 + 6 = 9
```

每个 H 的 STO-3G 是一个 s shell：

```text
H: [1s] -> 1 个 AO
```

两个 H：

```text
2 × 1 = 2
```

所以整个 H2O 的 AO 数是：

```text
O: 9
H: 2
total nao = 11
```

这里没有 d shell，所以 spherical 和 Cartesian 的 AO 数相同。

---

## 9. `ao_loc_nr()`：每个 shell 在 AO 序列中的起始位置

`nao_nr()` 只给出 AO 总数。但读积分矩阵、AO 切片时，还常需要知道：

```text
第 i 个 shell 对应 AO 编号的哪一段？
```

这就是 `ao_loc_nr()` 的作用。

源码里：

```python
def ao_loc_nr(mol, cart=None):
    if cart is None:
        cart = mol.cart
    if cart:
        return moleintor.make_loc(mol._bas, 'cart')
    else:
        return moleintor.make_loc(mol._bas, 'sph')
```

它返回一个长度为 `nbas + 1` 的数组。

如果：

```python
ao_loc = [0, 1, 4, 5]
```

就表示：

```text
shell 0 -> AO [0:1]
shell 1 -> AO [1:4]
shell 2 -> AO [4:5]
```

最后一个数就是 AO 总数：

```text
ao_loc[-1] == mol.nao_nr()
```

所以可以把 `ao_loc` 理解为 shell 到 AO 编号的索引表。

---

## 10. `nao_nr_range()`：某一段 shell 对应的 AO 范围

`nao_nr_range(mol, bas_id0, bas_id1)` 做的是：

```text
给定 shell 范围 [bas_id0, bas_id1)，返回对应 AO 范围 [nao_id0, nao_id1)
```

源码里：

```python
ao_loc = moleintor.make_loc(mol._bas[:bas_id1], 'sph')
nao_id0 = int(ao_loc[bas_id0])
nao_id1 = int(ao_loc[-1])
return nao_id0, nao_id1
```

例如文档里有：

```python
mol = gto.M(atom='O 0 0 0; C 0 0 1', basis='6-31g')
gto.nao_nr_range(mol, 2, 4)
# (2, 6)
```

意思是：

```text
第 2 到第 4 个 shell，对应 AO 编号 [2, 6)
```

这个函数在你想切某些 shell 对应的 AO block 时会很有用。

---

## 11. `nao_2c()`：spinor AO 数

`nao_2c()` 用于二分量 spinor GTO 的数量，和普通非相对论 AO 数不同。

源码里：

```python
def nao_2c(mol):
    l = mol._bas[:,ANG_OF]
    kappa = mol._bas[:,KAPPA_OF]
    dims = (l*4+2) * mol._bas[:,NCTR_OF]
    dims[kappa<0] = (l[kappa<0] * 2 + 2) * mol._bas[kappa<0,NCTR_OF]
    dims[kappa>0] = (l[kappa>0] * 2) * mol._bas[kappa>0,NCTR_OF]
    return int(dims.sum())
```

普通读 `mole.py` 第一阶段，可以先只记：

```text
nao_nr() -> 非相对论 AO 数
nao_2c() -> 二分量 spinor AO 数
```

`kappa` 是相对论/spinor basis 中的量子数标记；如果你现在主要读普通 Gaussian basis 和 SCF 主线，`nao_2c()` 可以先放后面。

---

## 12. `nao` 和矩阵维度的关系

`nao` 最直接的用处是决定 AO 表象下矩阵的维度。

例如：

```python
S = mol.intor('int1e_ovlp')
```

这个 overlap 矩阵的形状通常是：

```text
(nao, nao)
```

也就是：

$$
S_{\mu\nu} = \langle \chi_\mu | \chi_\nu \rangle
$$

其中 $\mu,\nu$ 都是在 AO basis 上跑的指标。

类似地：

```text
Hcore:      (nao, nao)
Density P:  (nao, nao)
Fock F:     (nao, nao)
MO coeff C: (nao, nmo)
```

所以 `nao` 不只是一个计数，它决定了后面 SCF 里很多矩阵的第一维和第二维。

---

## 13. 常见混淆点

### `nbas` 不是 `nao`

```text
nbas = shell 的数量
nao  = AO basis function 的数量
```

一个 shell 可以贡献多个 AO。

例如 p shell 是一个 shell，但通常贡献 3 个 AO。

### `nprim` 不是 `nao`

```text
nprim = primitive Gaussian 个数
nao   = contracted AO 个数
```

`nprim` 控制一个 contracted function 由多少个 primitive 线性组合而成，但不直接等于 AO 数。

### `nctr` 也不是完整的 `nao`

```text
nctr = contracted radial function 个数
nao  = nctr × 角向函数个数
```

例如一个 p shell 如果 `nctr=1`，仍然会贡献 3 个 AO。

### spherical 和 Cartesian 会影响 `nao`

从 d shell 开始：

```text
spherical d: 5 个 AO
Cartesian d: 6 个 AO
```

所以同一个 basis，在 `mol.cart=False` 和 `mol.cart=True` 下，`nao` 可能不同。

---

## 14. 读源码时的主线

读 `mole.py` 的 `nao` 相关函数，可以按这个顺序：

1. 先看 `_bas` 的 `ANG_OF / NPRIM_OF / NCTR_OF`
2. 再看 `nao_nr()`：用 `ANG_OF` 和 `NCTR_OF` 数 contracted spherical AO
3. 再看 `nao_cart()`：用 Cartesian 角向函数数
4. 再看 `npgto_nr()`：把 `NCTR_OF` 换成 `NPRIM_OF`，数 primitive GTO
5. 再看 `ao_loc_nr()`：从每个 shell 的 AO 数生成 AO 起始位置表
6. 最后再看 `nao_2c()`：相对论 spinor 情况

最值得记住的一句话是：

```text
nao_nr 数的是 contracted AO，不是 shell，不是 primitive。
```

再具体一点：

```text
每个 shell 的 AO 数 = 角向函数个数 × nctr
```