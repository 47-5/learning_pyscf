# 电子数、自旋、多重度：PySCF 中 `charge / spin / nelectron / nelec / multiplicity` 的关系

## 1. 先抓住三条核心公式

PySCF 里这组属性的核心关系是：

$$
N_e = \sum_A Z_A^{\mathrm{eff}} - q
$$

$$
\mathrm{spin} = N_\alpha - N_\beta = 2S
$$

$$
\mathrm{multiplicity} = 2S + 1 = \mathrm{spin} + 1
$$

其中：

- $N_e$ 是总电子数，对应 PySCF 的 `mol.nelectron`
- $q$ 是分子净电荷，对应 `mol.charge`
- $N_\alpha, N_\beta$ 对应 `mol.nelec` 返回的两个数
- $S$ 是总自旋量子数
- `mol.spin` 在 PySCF 里不是多重度，而是 $N_\alpha - N_\beta$

这和 Gaussian / ORCA 输入中的习惯不同。Gaussian / ORCA 通常写：

```text
charge multiplicity
```

而 PySCF 通常写：

```python
mol.charge = charge
mol.spin = multiplicity - 1
```

例如：

```text
Gaussian: 0 1  -> PySCF: charge=0, spin=0
Gaussian: 0 2  -> PySCF: charge=0, spin=1
Gaussian: 0 3  -> PySCF: charge=0, spin=2
Gaussian: 1 2  -> PySCF: charge=1, spin=1
```

---

## 2. `tot_electrons()`：总电子数从哪里来

`tot_electrons(mol)` 的主公式是：

$$
N_e = \sum_A Z_A^{\mathrm{eff}} - q
$$

在普通全电子计算中，$Z_A^{\mathrm{eff}}$ 就是原子核电荷。例如 H 是 1，O 是 8。

所以水分子中性体系：

$$
N_e = 8 + 1 + 1 - 0 = 10
$$

水分子阳离子：

$$
N_e = 8 + 1 + 1 - (+1) = 9
$$

源码中的逻辑大致是：

```python
if mol._atm.size != 0:
    nelectron = mol.atom_charges().sum()
elif mol._atom:
    nelectron = sum(charge(a[0]) for a in mol._atom)
else:
    nelectron = sum(charge(a[0]) for a in format_atom(mol.atom))
nelectron -= mol.charge
```

读法是：

- 如果 `_atm` 已经构建好，就从 `_atm` 里读原子电荷
- 如果 `_atm` 还没构建，但 `_atom` 已经有了，就从 `_atom` 读元素符号再查核电荷
- 如果连 `_atom` 也还没构建，就临时 `format_atom(mol.atom)` 再算
- 最后减去 `mol.charge`

这里要注意 ECP / pseudo。用了 ECP 或 pseudo 后，`_atm` 里的 `CHARGE_OF` 可能已经不是原始核电荷，而是有效核电荷或价电子相关的值。因此 `tot_electrons()` 在 `_atm` 构建后，会基于当前 `_atm` 看到的有效电荷来算。

---

## 3. `nelectron`：总电子数 property

`mol.nelectron` 是一个 property。它的 getter 逻辑是：

```python
@property
def nelectron(self):
    if self._nelectron is None:
        return self.tot_electrons()
    else:
        return self._nelectron
```

所以 `mol.nelectron` 有两种来源：

```text
用户没有显式设置 _nelectron
    -> 调用 tot_electrons() 根据原子和 charge 计算

用户显式设置过 mol.nelectron 或 mol.nelec
    -> 直接返回 self._nelectron
```

这就是为什么源码里有一个私有变量：

```python
self._nelectron = None
```

它的意思不是“电子数为 None”，而是：

```text
目前没有用户覆盖值，需要用 tot_electrons() 动态计算
```

---

## 4. `nelec`：返回 `(Nalpha, Nbeta)`，同时负责检查

`mol.nelec` 也是 property。它返回：

```python
(nalpha, nbeta)
```

数学关系是：

$$
N_\alpha + N_\beta = N_e
$$

$$
N_\alpha - N_\beta = \mathrm{spin}
$$

联立得到：

$$
N_\alpha = \frac{N_e + \mathrm{spin}}{2}
$$

$$
N_\beta = \frac{N_e - \mathrm{spin}}{2}
$$

PySCF 代码里写成：

```python
ne = self.nelectron
nalpha = (ne + self.spin) // 2
nbeta = nalpha - self.spin
```

这里的 `//` 是整数除法。真正的合法性检查在后面：

```python
if nalpha + nbeta != ne:
    raise RuntimeError(...)
```

这一步检查的是：

```text
nelectron 和 spin 的奇偶性是否匹配
```

例如中性水分子：

```text
ne = 10
spin = 0
nalpha = (10 + 0) / 2 = 5
nbeta = 5
合法
```

如果你错误地设成：

```text
ne = 10
spin = 1
```

理论上会得到：

$$
N_\alpha = \frac{10+1}{2}=5.5
$$

$$
N_\beta = \frac{10-1}{2}=4.5
$$

半个电子不合法。代码里因为用了 `//`，会先得到整数近似，然后用 `nalpha + nbeta != ne` 把这个错误抓出来。

所以 `nelec` 不只是“返回电子数拆分”，它还负责：

```text
检查 spin 和总电子数是否自洽
```

---

## 5. 为什么 `build()` 里会写一行孤零零的 `self.nelec`

`build()` 里有这样一段：

```python
if self.spin is None:
    self.spin = self.nelectron % 2
else:
    self.nelec
```

这行 `self.nelec` 看起来没有用返回值，但它不是无意义代码。

因为 `nelec` 是 `@property`，访问它就会执行 getter：

```python
@property
def nelec(self):
    ...
```

也就是说：

```python
self.nelec
```

等价于：

```text
请运行 nelec getter，并在里面检查 nelectron 和 spin 是否自洽
```

如果自洽，什么事都不发生；如果不自洽，就在这里抛出错误。

所以这行代码的目的不是拿返回值，而是触发校验。

---

## 6. `spin is None`：让 PySCF 自动猜自旋

如果用户设置：

```python
mol.spin = None
```

或者：

```python
mol.multiplicity = None
```

那么 `build()` 会进入：

```python
self.spin = self.nelectron % 2
```

这表示：

```text
偶数电子 -> spin = 0
奇数电子 -> spin = 1
```

也就是给出最小的、和电子数奇偶性匹配的自旋设定。

例如：

```text
H2O: ne = 10 -> spin = 0 -> singlet
H atom: ne = 1 -> spin = 1 -> doublet
H2O+: ne = 9 -> spin = 1 -> doublet
```

这只是一个默认猜测，不等于程序知道真实基态自旋。比如 O2 的真实基态是 triplet，你通常应该显式写：

```python
mol.spin = 2
```

---

## 7. `multiplicity`：给熟悉 Gaussian / ORCA 的人用的桥

PySCF 里有：

```python
@property
def multiplicity(self):
    return self.spin + 1

@multiplicity.setter
def multiplicity(self, x):
    if x is None:
        self.spin = None
    else:
        self.spin = x - 1
```

所以你可以用多重度来设置自旋：

```python
mol.multiplicity = 1  # spin = 0
mol.multiplicity = 2  # spin = 1
mol.multiplicity = 3  # spin = 2
```

但在很多 PySCF 示例里，更常直接写：

```python
mol.spin = 2
```

这就是 triplet，而不是 doublet。

最容易记错的是：

```text
PySCF spin = multiplicity - 1
Gaussian/ORCA 第二个数字 = multiplicity
```

---

## 8. `nelec` setter：反过来用 alpha/beta 设置体系

`nelec` 也有 setter：

```python
@nelec.setter
def nelec(self, neleca_nelecb):
    neleca, nelecb = neleca_nelecb
    self._nelectron = neleca + nelecb
    self.spin = neleca - nelecb
```

所以如果你写：

```python
mol.nelec = (5, 5)
```

那么 PySCF 会设置：

```text
_nelectron = 10
spin = 0
```

如果你写：

```python
mol.nelec = (9, 7)
```

那么：

```text
_nelectron = 16
spin = 2
multiplicity = 3
```

这对需要手动控制 alpha/beta 电子数的情况很方便。

---

## 9. `ms`：自旋量子数 S

PySCF 还有一个 property：

```python
@property
def ms(self):
    if self.spin % 2 == 0:
        return self.spin // 2
    else:
        return self.spin * .5
```

这里的 `ms` 可以理解为：

$$
S = \frac{\mathrm{spin}}{2}
$$

因此：

```text
spin = 0 -> ms = 0
spin = 1 -> ms = 0.5
spin = 2 -> ms = 1
```

对应多重度：

$$
\mathrm{multiplicity} = 2S + 1
$$

---

## 10. 几个具体例子

| 体系 | charge | nelectron | PySCF spin | multiplicity | nelec `(Nalpha, Nbeta)` |
|---|---:|---:|---:|---:|---:|
| H2O singlet | 0 | 10 | 0 | 1 | `(5, 5)` |
| H atom doublet | 0 | 1 | 1 | 2 | `(1, 0)` |
| H2O+ doublet | +1 | 9 | 1 | 2 | `(5, 4)` |
| O2 triplet | 0 | 16 | 2 | 3 | `(9, 7)` |
| 闭壳层阴离子例子 | -1 | 奇偶取决于体系 | 需要手动检查 | 取决于 spin | 取决于 nelectron/spin |

一个错误例子：

```python
mol.atom = 'O 0 0 0; H 0 0 1; H 0 1 0'
mol.charge = 0
mol.spin = 1
mol.build()
```

中性水分子有 10 个电子，但 `spin=1` 要求 alpha/beta 相差 1。这样会导致半个电子，所以 `build()` 中访问 `self.nelec` 时会抛错。

---

## 11. 这组属性的调用图

可以把调用关系记成：

```text
mol.charge
    -> 影响 tot_electrons()
    -> 影响 mol.nelectron

mol.nelectron getter
    -> 如果 _nelectron is None: 调用 tot_electrons()
    -> 否则返回 _nelectron

mol.nelec getter
    -> 读取 mol.nelectron
    -> 读取 mol.spin
    -> 计算 (Nalpha, Nbeta)
    -> 检查 electron number 和 spin 是否自洽

mol.multiplicity setter
    -> 设置 mol.spin = multiplicity - 1

mol.nelec setter
    -> 设置 mol._nelectron = Nalpha + Nbeta
    -> 设置 mol.spin = Nalpha - Nbeta

build()
    -> 如果 spin is None: spin = nelectron % 2
    -> 否则访问 self.nelec，触发一致性检查
```

---

## 12. 读源码时的主线

读 `mole.py` 这块时，不要把这些 property 看成孤立函数。它们是一组互相配合的接口：

```text
charge 表示分子净电荷，负责改总电子数
nelectron 负责给出总电子数
spin 负责指定 alpha/beta 差值
nelec 负责把总电子数和 spin 拆成 alpha/beta，并检查一致性
multiplicity 负责把化学多重度翻译成 PySCF spin
build 负责在真正构建 Mole 时做最后检查
```

最值得记住的一句话是：

```text
PySCF 的 spin = Nalpha - Nbeta = 2S，不是 multiplicity。
```

以及：

```text
build() 里单独访问 self.nelec，是为了触发 property getter 中的自洽性检查。
```