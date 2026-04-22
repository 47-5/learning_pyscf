# `make_env` 读法：用 H2O 例子理解 `_atm`、`_bas`、`_env`

## 1. 先说结论

`make_env` 的作用不是做化学计算，而是把 PySCF 里比较“人类可读”的原子和基组信息，打包成 `libcint` 可以直接使用的三块数组：

```text
_atm  = 原子目录
_bas  = 壳层目录
_env  = 真正的浮点数据仓库
```

可以把它理解成：

- `_atm` 说“原子是谁、坐标放在哪”
- `_bas` 说“壳层是什么、挂在哪个原子上、参数放在哪”
- `_env` 真正存坐标、指数、收缩系数等浮点数

## 2. 本例设置：一个简化的 H2O

为了让 `_atm/_bas/_env` 里的指针数字尽量好看，这里不用真实实验几何，而是直接取一个简单的 Bohr 坐标例子：

```python
atom = [
    ['O', [0.0,  0.0, 0.0]],
    ['H', [0.0, -1.5, 1.0]],
    ['H', [0.0,  1.5, 1.0]],
]

basis = {
    'O': '6-31g',
    'H': 'sto-3g',
}
```

这里假设：

- 坐标已经是 `format_atom` 处理后的内部格式
- 单位已经是 Bohr
- `pre_env` 只有默认的前 20 个保留槽位
- `basis` 字典按 `O` 然后 `H` 的顺序进入 `make_env`

## 3. H 和 O 的基组先长什么样

### H: STO-3G

H 的 STO-3G 在基组文件里对应一个 s shell：

```text
H    S
  3.42525091   0.15432897
  0.62391373   0.53532814
  0.16885540   0.44463454
```

内部上可以看成：

```python
[0,
 [3.42525091, 0.15432897],
 [0.62391373, 0.53532814],
 [0.16885540, 0.44463454]]
```

所以：

- `l = 0`
- `nprim = 3`
- `nctr = 1`
- 每个 H 原子只会贡献 1 个 s shell

### O: 6-31G

O 的 6-31G 在文件里对应三段：

```text
O    S
  5484.6717000   0.0018311
   825.2349500   0.0139501
   188.0469600   0.0684451
    52.9645000   0.2327143
    16.8975700   0.4701930
     5.7996353   0.3585209
O    SP
    15.5396160  -0.1107775   0.0708743
     3.5999336  -0.1480263   0.3397528
     1.0137618   1.1307670   0.7271586
O    SP
     0.2700058   1.0000000   1.0000000
```

这里最重要的是：`SP` 会被拆成一个 s shell 和一个 p shell。

所以 O 的 6-31G 最终会变成 5 个 shell：

1. 一个 6 primitive 的 s shell
2. 一个 3 primitive 的 s shell
3. 一个 3 primitive 的 p shell
4. 一个 1 primitive 的 s shell
5. 一个 1 primitive 的 p shell

也就是常说的：

```text
(10s,4p) -> [3s,2p]
```

## 4. `pre_env` 从哪里开始

PySCF 在 `_env` 开头先预留了 20 个控制槽位：

```python
PTR_ENV_START = 20
```

所以如果没有别的额外数据，`make_env` 一开始相当于：

```python
pre_env = [0.0] * 20
ptr_env = 20
```

这意味着：

- `_env[0:20]` 先留给 libcint 控制参数
- 分子自己的坐标、指数、系数从 `_env[20]` 开始往后写

## 5. 第一步：先处理三个原子，得到 `_atm`

`make_env` 的前半段会先循环 `atoms`，每个原子调用一次 `make_atm_env`。

### O 原子

O 在 `(0,0,0)`，而当前 `ptr_env = 20`。

于是 `make_atm_env` 生成：

```python
env0 = [0.0, 0.0, 0.0, 0.0] # [x, y, z, zeta(核当成点电荷时为0)]
atm0 = [8, 20, 1, 23, 0, 0]
```

解释：

```text
[CHARGE_OF, PTR_COORD, NUC_MOD_OF, PTR_ZETA, PTR_FRAC_CHARGE, PTR_RADIUS]
=
[8,         20,        1,          23,       0,               0]
```

也就是：

- 核电荷是 8
- 坐标在 `_env[20:23]`
- 核模型是点核 `NUC_POINT = 1`
- zeta 在 `_env[23]`

写完 O 以后，`ptr_env = 24`。

### 第一个 H 原子

第一个 H 在 `(0,-1.5,1.0)`，此时 `ptr_env = 24`。

得到：

```python
env1 = [0.0, -1.5, 1.0, 0.0]
atm1 = [1, 24, 1, 27, 0, 0]
```

写完以后，`ptr_env = 28`。

### 第二个 H 原子

第二个 H 在 `(0,1.5,1.0)`，此时 `ptr_env = 28`。

得到：

```python
env2 = [0.0, 1.5, 1.0, 0.0]
atm2 = [1, 28, 1, 31, 0, 0]
```

写完以后，`ptr_env = 32`。

### 这时 `_atm` 是什么

所以原子部分结束后：

```python
_atm = [
    [8, 20, 1, 23, 0, 0],
    [1, 24, 1, 27, 0, 0],
    [1, 28, 1, 31, 0, 0],
]
```

而 `_env` 已经占用了：

```text
0:19   预留控制槽位
20:23  O 的 [x, y, z, zeta]
24:27  H1 的 [x, y, z, zeta]
28:31  H2 的 [x, y, z, zeta]
```

## 6. 第二步：按元素构造 basis 模板

接下来 `make_env` 会按元素循环 `basis.items()`，调用 `make_bas_env`。

注意：

**这里不是给每个原子都存一遍基组数字，而是同一种元素先存一份模板。**

### O 的 6-31G 模板

此时 `ptr_env = 32`，开始写 O 的基组数据。

O 一共有 5 个 shell，因此 `_basdic['O']` 会先得到 5 行模板。

#### O shell 1: 6 primitive 的 s shell

```text
l = 0, nprim = 6, nctr = 1
ptr_exp   = 32
ptr_coeff = 38
env 长度  = 6 + 6 = 12
```

对应 `_bas` 模板行：

```python
[0, 0, 6, 1, 0, 32, 38, 0]
```

写完后，下一个空位是 `44`。

#### O shell 2: 由第一段 `SP` 拆出来的 s shell

```python
[0, 0, 3, 1, 0, 44, 47, 0]
```

长度 `3 + 3 = 6`，写完后下一个空位是 `50`。

#### O shell 3: 由第一段 `SP` 拆出来的 p shell

```python
[0, 1, 3, 1, 0, 50, 53, 0]
```

长度 `3 + 3 = 6`，写完后下一个空位是 `56`。

#### O shell 4: 由第二段 `SP` 拆出来的 s shell

```python
[0, 0, 1, 1, 0, 56, 57, 0]
```

长度 `1 + 1 = 2`，写完后下一个空位是 `58`。

#### O shell 5: 由第二段 `SP` 拆出来的 p shell

```python
[0, 1, 1, 1, 0, 58, 59, 0]
```

长度 `1 + 1 = 2`，写完后下一个空位是 `60`。

所以 O 这套模板一共占用了 `_env[32:60]`，总长度 28。

### H 的 STO-3G 模板

现在 `ptr_env = 60`。

H 只有 1 个 s shell：

```python
[0, 0, 3, 1, 0, 60, 63, 0]
```

对应：

- `ptr_exp = 60`
- `ptr_coeff = 63`
- 长度 `3 + 3 = 6`

写完后，下一个空位是 `66`。

所以 H 的模板占用 `_env[60:66]`。

## 7. 第三步：把模板分配给具体原子

此时 `_basdic` 里已经有：

- `_basdic['O']`：5 行 O 的壳层模板
- `_basdic['H']`：1 行 H 的壳层模板

接着 `make_env` 再循环一遍原子，把模板分配给每个具体原子，并把 `ATOM_OF` 改成真实原子编号。

### O 原子（atom 0）

O 拿到 5 行模板，并把 `ATOM_OF` 改成 0：

```python
[0, 0, 6, 1, 0, 32, 38, 0]
[0, 0, 3, 1, 0, 44, 47, 0]
[0, 1, 3, 1, 0, 50, 53, 0]
[0, 0, 1, 1, 0, 56, 57, 0]
[0, 1, 1, 1, 0, 58, 59, 0]
```

### 第一个 H 原子（atom 1）

H 拿到 1 行模板，并把 `ATOM_OF` 改成 1：

```python
[1, 0, 3, 1, 0, 60, 63, 0]
```

### 第二个 H 原子（atom 2）

再拿到同一行模板的拷贝，并把 `ATOM_OF` 改成 2：

```python
[2, 0, 3, 1, 0, 60, 63, 0]
```

注意最后这两行非常关键：

- 第一个 H 和第二个 H 的 `ATOM_OF` 不同
- 但 `PTR_EXP = 60`、`PTR_COEFF = 63` 相同

这说明：

**两个 H 共用同一份 STO-3G 的 exponent / coefficient 数据，只是壳层中心不同。**

## 8. 最终 `_bas` 长什么样

所以最终 `_bas` 一共 7 行：

```python
_bas = [
    [0, 0, 6, 1, 0, 32, 38, 0],
    [0, 0, 3, 1, 0, 44, 47, 0],
    [0, 1, 3, 1, 0, 50, 53, 0],
    [0, 0, 1, 1, 0, 56, 57, 0],
    [0, 1, 1, 1, 0, 58, 59, 0],
    [1, 0, 3, 1, 0, 60, 63, 0],
    [2, 0, 3, 1, 0, 60, 63, 0],
]
```

逐行解释就是：

- 前 5 行是 O 的 5 个 shell
- 第 6 行是第一个 H 的 1s 壳层
- 第 7 行是第二个 H 的 1s 壳层

## 9. 最终 `_env` 的布局怎么读

在这个例子里，你不一定要把所有归一化后的 coefficient 小数完全背下来，更重要的是看懂布局。

因此 `_env` 可以先读成：

```text
0:19    预留控制槽位
20:23   O 的 [x, y, z, zeta]
24:27   H1 的 [x, y, z, zeta]
28:31   H2 的 [x, y, z, zeta]
32:43   O shell 1 的 exponent + normalized coeff
44:49   O shell 2 的 exponent + normalized coeff
50:55   O shell 3 的 exponent + normalized coeff
56:57   O shell 4 的 exponent + normalized coeff
58:59   O shell 5 的 exponent + normalized coeff
60:65   H shell 的 exponent + normalized coeff
```

注意：

- `make_bas_env` 里会先做 primitive 归一化
- 如果 `NORMALIZE_GTO = True`，还会继续做 contracted AO 归一化

所以 `_env` 里的系数不是用户在基组文件里看到的原始系数，而是 `libcint` 真正要用的归一化后系数。

## 10. `libcint` 如何利用这些数组

举两个例子最直观。

### 例 1：读取第二个 H 的壳层

看 `_bas` 最后一行：

```python
[2, 0, 3, 1, 0, 60, 63, 0]
```

可读成：

- `ATOM_OF = 2`：这个壳层挂在第 2 号原子，也就是第二个 H
- `ANG_OF = 0`：s shell
- `NPRIM_OF = 3`
- `NCTR_OF = 1`
- `PTR_EXP = 60`
- `PTR_COEFF = 63`

再看 `_atm[2]`：

```python
[1, 28, 1, 31, 0, 0]
```

可读成：

- 核电荷是 1
- 坐标在 `_env[28:31]`

于是 `libcint` 就能拼出这个壳层的全部必要信息：

- 中心坐标：`_env[28:31]`
- exponents：`_env[60:63]`
- coefficients：`_env[63:66]`

### 例 2：读取 O 的第一个 p shell

看 `_bas` 第 3 行：

```python
[0, 1, 3, 1, 0, 50, 53, 0]
```

可读成：

- 这个壳层属于 O
- 角动量 `l = 1`，所以是 p shell
- 有 3 个 primitive
- 有 1 个 contraction
- exponents 在 `_env[50:53]`
- coefficients 在 `_env[53:56]`

再结合 `_atm[0]` 里的 O 坐标，积分库就知道这个 p shell 以 O 为中心展开。

## 11. 这个例子里 AO 数有多少

如果按非相对论 spherical AO 来看：

- O 的 3 个 s shell 贡献 3 个 AO
- O 的 2 个 p shell 各贡献 3 个 AO，共 6 个 AO
- 两个 H 各 1 个 s AO，共 2 个 AO

所以总 AO 数是：

```text
3 + 6 + 2 = 11
```

这和你对 H2O / 6-31G + STO-3G 这种组合的直觉是一致的。

## 12. 这一页真正想记住什么

这个 H2O 例子里，最值得记住的不是每个数字，而是下面三句话：

```text
_atm 记录“原子是谁、坐标在哪里”
_bas 记录“壳层是什么、挂在哪、参数在哪里”
_env 真正存“坐标、指数、系数”等浮点数据
```

再进一步就是：

```text
同一种元素的基组参数通常只存一份模板
每个具体原子只是复制壳层记录，并修改 ATOM_OF
```

这正是 `make_env` 的核心设计。

## 13. MYNOTE：下一步该怎么看

看完这一页之后，最适合继续读的是：

1. `make_atm_env`
2. `make_bas_env`
3. `make_env`
4. `atom_coord()`
5. `bas_exp()`
6. `bas_ctr_coeff()`

因为前面三个函数负责“打包进去”，后面三个函数负责“再从 `_atm/_bas/_env` 里读出来”。

如果这两头都看通了，`mole.py` 里关于基组与原子中心的底层表示就真的连起来了。

## 14. 附注：libcint 输入数组的布局说明

``````python
# libcint 输入数组的布局说明。
# _atm[atom_id]：每个原子一行，存原子的整数元信息。
# _bas[shell_id]：每个壳层一行，存壳层的整数元信息。
# _env：真正存放坐标、指数、系数等浮点数据的大数组。
# *_OF 表示 _atm/_bas 的列号；PTR_* 表示指向 _env 的位置。
# 注意：有些槽位会在 ECP / spinor 路径中复用。
CHARGE_OF  = 0  # _atm 列：核电荷 Z
PTR_COORD  = 1  # _atm 列：该原子 xyz 坐标在 _env 中的起始位置
NUC_MOD_OF = 2  # _atm 列：核模型标记（点核 / Gaussian 核 / 分数电荷 / ECP）
PTR_ZETA   = 3  # _atm 列：核模型 zeta 在 _env 中的位置
PTR_FRAC_CHARGE = 4  # _atm 列：MM/分数电荷在 _env 中的位置
PTR_RADIUS = 5  # _atm 列：预留的半径指针
ATM_SLOTS  = 6  # _atm 每一行的槽位数

ATOM_OF    = 0  # _bas 列：该壳层属于哪个原子
ANG_OF     = 1  # _bas 列：角动量 l
NPRIM_OF   = 2  # _bas 列：primitive Gaussian 的个数
NCTR_OF    = 3  # _bas 列：contracted function 的个数
RADI_POWER = 3  # _ecpbas 列：径向幂次槽位
KAPPA_OF   = 4  # _bas 列：spinor/相对论壳层使用的 kappa 标记
SO_TYPE_OF = 4  # _ecpbas 列：自旋轨道 ECP 类型槽位
PTR_EXP    = 5  # _bas 列：指数 exponent 在 _env 中的起始位置
PTR_COEFF  = 6  # _bas 列：收缩系数 coefficient 在 _env 中的起始位置
BAS_SLOTS  = 8  # _bas 每一行的槽位数

# _env 开头预留的控制槽位。
PTR_EXPCUTOFF   = 0   # _env 槽位：libcint 的指数截断阈值
PTR_COMMON_ORIG = 1   # _env[1:4]：诸如 r 这类算符共用的原点
PTR_RINV_ORIG   = 4   # _env[4:7]：1/|r-R| 算符使用的原点
PTR_RINV_ZETA   = 7   # _env 槽位：rinv 算符中 Gaussian 核模型的 zeta
PTR_RANGE_OMEGA = 8   # _env 槽位：range-separated Coulomb 的 omega
PTR_F12_ZETA    = 9   # _env 槽位：F12 geminal 指数
PTR_GTG_ZETA    = 10  # _env 槽位：GTG geminal 指数
NGRIDS          = 11  # _env 槽位：int1e_grids 中的网格点个数
PTR_GRIDS       = 12  # _env 槽位：展开后网格坐标在 _env 中的位置
AS_RINV_ORIG_ATOM = 17  # _env 槽位：rinv/ECP 导数辅助用到的原子编号
AS_ECPBAS_OFFSET = 18   # _env 槽位：_ecpbas 相对 _bas 的起始偏移
AS_NECPBAS      = 19    # _env 槽位：ECP basis 记录条数
PTR_ENV_START   = 20    # 分子自身数据在 _env 中真正开始写入的位置

# 存在 _atm[:, NUC_MOD_OF] 里的核模型枚举值。
NUC_POINT = 1  # 点核模型
NUC_GAUSS = 2  # 有限大小的 Gaussian 核模型
NUC_FRAC_CHARGE = 3  # MM / 分数电荷粒子
NUC_ECP = 4  # 使用 ECP / pseudo 描述的原子
``````

