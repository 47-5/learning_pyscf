# GTO 归一化：primitive、contracted 与 PySCF 中的处理

## 1. 这份笔记想解决什么问题

在读 `pyscf/gto/mole.py` 时，归一化相关的代码主要集中在这几行：

```python
cs = numpy.einsum('pi,p->pi', cs, gto_norm(angl, es))
if NORMALIZE_GTO:
    cs = _nomalize_contracted_ao(angl, es, cs)
```

如果不先把“primitive 归一化”和“contracted 归一化”分开，很容易把下面几件事混在一起：

1. 单个 primitive Gaussian 自己要不要归一化
2. 多个 primitive 线性组合成一个 contracted function 之后要不要再归一化
3. PySCF 到底是把归一化因子放在 basis function 里，还是吸收到 coefficient 里
4. `bas_ctr_coeff()` 返回的系数，和 basis 文件中的原始数字是不是一回事

这份笔记主要讨论 `mole.py` 里 `gto_norm()` 和 `_nomalize_contracted_ao()` 对应的那部分内容。

---

## 2. primitive GTO 的归一化

先以球谐视角写一个中心在原子 $A$ 上的 primitive GTO：

$$
\phi_{plm}^{A}(\mathbf r)
=
N_l(\alpha_p)
\,r_A^l e^{-\alpha_p r_A^2}
\,Y_{lm}(\hat{\mathbf r}_A)
$$

其中：

- $p$ 表示第几个 primitive exponent
- $l,m$ 是角动量量子数
- $r_A = |\mathbf r - \mathbf R_A|$
- $Y_{lm}$ 是角向部分
- $N_l(\alpha_p)$ 是我们要找的归一化因子

这里最重要的一点是：在这种写法下，角向部分通常取为已经归一化的球谐函数，即

$$
\int |Y_{lm}(\hat{\mathbf r})|^2\,d\Omega = 1
$$

所以 primitive GTO 的归一化条件主要落在径向部分上：

$$
\langle \phi_{plm}^{A} | \phi_{plm}^{A} \rangle = 1
$$

代入后得到：

$$
1
=
N_l(\alpha_p)^2
\int_0^{\infty}
\left(r_A^l e^{-\alpha_p r_A^2}\right)^2
r_A^2\,dr_A
$$

也就是

$$
1
=
N_l(\alpha_p)^2
\int_0^{\infty} r_A^{2l+2} e^{-2\alpha_p r_A^2}\,dr_A
$$

因此

$$
N_l(\alpha_p)
=
\left(
\int_0^{\infty} r^{2l+2} e^{-2\alpha_p r^2}\,dr
\right)^{-1/2}
$$

PySCF 在源码里把这个积分包装成了：

$$
\mathrm{gaussian\_int}(n,\beta)
=
\int_0^{\infty} x^n e^{-\beta x^2}\,dx
$$

于是 primitive 归一化因子就可以写成：

$$
N_l(\alpha)
=
\frac{1}{\sqrt{\mathrm{gaussian\_int}(2l+2, 2\alpha)}}
$$

这正是 `gto_norm(l, expnt)` 在做的事情。

PySCF 源码中的对应关系是：

```python
def gaussian_int(n, alpha):
    return scipy.special.gamma((n + 1) * .5) / (2. * alpha**((n + 1) * .5))

def gto_norm(l, expnt):
    return 1/numpy.sqrt(gaussian_int(l*2+2, 2*expnt))
```

也就是说，`gto_norm()` 返回的是“把径向 primitive Gaussian 归一化到 1 所需要乘上的因子”。

---

## 3. 一个重要提醒：这里主要是径向归一化

在 `make_bas_env()` 这条路径里，`gto_norm(l, expnt)` 只显式依赖于 $l$ 和 $\alpha$，并不区分 $m$，也不区分某个 Cartesian 分量到底是 $x^2y$ 还是 $xyz$。

所以从阅读 `mole.py` 的角度，你可以先把它理解为：

- 这里优先处理的是 **radial part** 的归一化
- 角向部分的规范化约定由球谐/Cartesian 约定和 `libcint` 进一步处理

对于你当前在读的 `make_bas_env()`，先把注意力集中在“primitive 的径向归一化”和“contracted function 的整体归一化”上就够了。

---

## 4. 两种等价的记账方式

这是理解 PySCF 实现时最关键的一步。

设未归一化的 primitive 径向函数记为：

$$
g_p(r) = r^l e^{-\alpha_p r^2}
$$

则归一化后的 primitive 可以写成：

$$
\tilde g_p(r) = N_p g_p(r)
$$

其中 $N_p = N_l(\alpha_p)$。

现在考虑一个 contracted radial function：

$$
R_k(r) = \sum_{p=1}^{n_{\mathrm{prim}}} c_{pk}\,\tilde g_p(r)
$$

这里有两种完全等价的写法。

### 写法 A：把归一化因子放在 basis function 里

$$
R_k(r)
=
\sum_p c_{pk}\,\tilde g_p(r)
=
\sum_p c_{pk}\,N_p g_p(r)
$$

这里：

- basis function 是 $\tilde g_p$
- coefficient 是原始的 $c_{pk}$

### 写法 B：把归一化因子吸收到 coefficient 里

定义新的系数：

$$
d_{pk} = N_p c_{pk}
$$

那么就可以改写成：

$$
R_k(r) = \sum_p d_{pk}\,g_p(r)
$$

这里：

- basis function 是未归一化的 $g_p$
- coefficient 变成了已经吸收 primitive 归一化因子的 $d_{pk}$

这两种写法在数学上是等价的。

**PySCF 在 `make_bas_env()` 里采用的是写法 B。**

也就是说，PySCF 会先把每个 primitive 的归一化因子乘到系数矩阵里，而不是把它保留成一个独立的“basis function 前因子”。

---

## 5. 为什么 contracted function 还要再归一化一次

就算每个 primitive 自己都已经归一化了，多个 primitive 的线性组合也不一定自动归一化。

原因很简单：不同 primitive 彼此通常 **不正交**。

如果只有一个 contraction，写成：

$$
R(r) = c_1 \tilde g_1(r) + c_2 \tilde g_2(r)
$$

那么它的范数平方是：

$$
\langle R|R\rangle
=
c_1^2 \langle \tilde g_1|\tilde g_1\rangle
+ c_2^2 \langle \tilde g_2|\tilde g_2\rangle
+ 2 c_1 c_2 \langle \tilde g_1|\tilde g_2\rangle
$$

因为

$$
\langle \tilde g_1|\tilde g_2\rangle \neq 0
$$

所以一般不能靠

$$
c_1^2 + c_2^2 = 1
$$

来完成归一化。

也就是说，**contracted function 的归一化不是普通欧氏长度问题，而是带重叠矩阵的二次型问题。**

---

## 6. overlap 矩阵是怎么来的

如果采用上面的写法 B，也就是把 primitive 归一化因子已经吸收到系数中，那么基函数仍然写成未归一化的

$$
g_p(r) = r^l e^{-\alpha_p r^2}
$$

此时 primitive 间的径向 overlap 矩阵元素就是：

$$
S_{pq}^{\mathrm{raw}}
=
\int_0^{\infty} g_p(r) g_q(r) r^2\,dr
=
\int_0^{\infty} r^{2l+2} e^{-(\alpha_p+\alpha_q)r^2}\,dr
$$

所以

$$
S_{pq}^{\mathrm{raw}} = \mathrm{gaussian\_int}(2l+2, \alpha_p+\alpha_q)
$$

这就是 `_nomalize_contracted_ao()` 里：

```python
ee = es.reshape(-1,1) + es.reshape(1,-1)
ee = gaussian_int(l*2+2, ee)
```

的数学含义。注：es: exponent 向量，长度为 nprim

这里：

- $es[p] = \alpha_p$
- $ee[p,q] = \alpha_p + \alpha_q$
- $\mathrm{gaussian\_int}(2l+2, ee)[p,q] = S_{pq}^{\mathrm{raw}}$

如果想和代码逐字对应，也可以把第三条读成：代码里写作 `gaussian_int(l*2+2, ee)`，数学上就是 $\mathrm{gaussian\_int}(2l+2, ee)$。

因此 $ee$ 本质上就是 primitive 径向函数之间的 overlap 矩阵。

---

## 7. contracted 归一化的矩阵写法

把 primitive 归一化已经吸收到系数之后，记某个 contraction 的系数列向量为：

$$
\mathbf d_k
=
\begin{bmatrix}
d_{1k} \\
d_{2k} \\
\vdots \\
d_{n_{\mathrm{prim}}k}
\end{bmatrix}
$$

对应的 contracted radial function 是：

$$
R_k(r) = \sum_p d_{pk} g_p(r)
$$

那么它的范数平方就是：

$$
\langle R_k | R_k \rangle
=
\mathbf d_k^T \mathbf S^{\mathrm{raw}} \mathbf d_k
$$

如果你暂时不熟悉这种矩阵写法，可以先把

$$
\mathbf d_k^T \mathbf S^{\mathrm{raw}} \mathbf d_k
$$

理解成：**在非正交 primitive 基底里计算“长度平方”的正确方式。**

它的双重求和展开、积分形式推导，以及两 primitive 特例，我都放到了文末“附录 A”里，正文这里先只保留主线。

因此，使它归一化所需的缩放因子为：

$$
s_k
=
\frac{1}{\sqrt{\mathbf d_k^T \mathbf S^{\mathrm{raw}} \mathbf d_k}}
$$

归一化后的系数列向量为：

$$
\mathbf d_k^{\mathrm{norm}} = s_k \mathbf d_k
$$

如果一个 shell 有多个 contraction，那么把所有列拼成矩阵：

$$
\mathbf D
=
\begin{bmatrix}
d_{11} & d_{12} & \cdots & d_{1n_{\mathrm{ctr}}} \\
d_{21} & d_{22} & \cdots & d_{2n_{\mathrm{ctr}}} \\
\vdots & \vdots & \ddots & \vdots \\
d_{n_{\mathrm{prim}}1} & d_{n_{\mathrm{prim}}2} & \cdots & d_{n_{\mathrm{prim}}n_{\mathrm{ctr}}}
\end{bmatrix}
\in \mathbb R^{n_{\mathrm{prim}}\times n_{\mathrm{ctr}}}
$$

那么每一列都要各自乘上自己的归一化因子：

$$
\mathbf D_{\mathrm{norm}}
=
\mathbf D
\begin{bmatrix}
s_1 & 0 & \cdots & 0 \\
0 & s_2 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & s_{n_{\mathrm{ctr}}}
\end{bmatrix}
$$

这就是 `_nomalize_contracted_ao()` 最后这句的意义：

```python
return numpy.einsum('pi,i->pi', cs, s1)
```

也就是：

- `cs` 的每一列是一组 contraction
- `s1[i]` 是第 `i` 列的整体归一化因子
- 把第 `i` 列整体乘上 `s1[i]`

---

## 8. PySCF 代码到底在做什么

现在回头看 `make_bas_env()` 里的关键三步：

### 第一步：读出 exponent 和 coefficient

```python
es = b_coeff[:,0]
cs = b_coeff[:,1:]
```

这里：

- `es.shape = (nprim,)`
- `cs.shape = (nprim, nctr)`

在这一刻，`cs` 还是“从 basis 数据读来的原始 contraction coefficient 矩阵”。

### 第二步：把 primitive 归一化因子吸收到系数里

```python
cs = numpy.einsum('pi,p->pi', cs, gto_norm(angl, es))
```

若记

$$
\mathbf N = \mathrm{diag}(N_1, N_2, \cdots, N_{n_{\mathrm{prim}}})
$$

那么这一步就是：

$$
\mathbf C_{\mathrm{raw}}
\longrightarrow
\mathbf D = \mathbf N \, \mathbf C_{\mathrm{raw}}
$$

也就是把每一行分别乘上对应 primitive 的归一化因子。（行对应的一个源GTO）

### 第三步：对每一列 contraction 再做整体归一化

```python
ee = gaussian_int(l*2+2, es.reshape(-1,1) + es.reshape(1,-1))
s1 = 1. / numpy.sqrt(numpy.einsum('pi,pq,qi->i', cs, ee, cs))
cs = numpy.einsum('pi,i->pi', cs, s1)
```

若此时 `cs` 已经是上面的 $\mathbf D$，那么：

$$
s_k = \frac{1}{\sqrt{\mathbf D_{:k}^T \mathbf S^{\mathrm{raw}} \mathbf D_{:k}}}
$$

最后得到：

$$
\mathbf C_{\mathrm{env}}
=
\mathbf N \, \mathbf C_{\mathrm{raw}} \, \mathbf S_c
$$

其中

$$
\mathbf S_c = \mathrm{diag}(s_1, s_2, \cdots, s_{n_{\mathrm{ctr}}})
$$

也就是说，`_env` 里最终存下来的系数矩阵，已经同时包含了：

1. primitive 归一化
2. contracted 归一化

因此它是 **libcint 直接可用** 的系数矩阵。

*注意第二步和第三步einsum的写法不同，一个是左乘（对行），一个是右乘（对列）*

---

## 9. 一个非常容易混淆的点：归一化不等于正交化

即使做完了 contracted normalization，也只是保证每一列满足：

$$
\langle R_k | R_k \rangle = 1
$$

并不保证不同列之间满足：

$$
\langle R_k | R_{k'} \rangle = 0
\qquad (k \neq k')
$$

换成矩阵语言，就是：

$$
\mathbf D_{\mathrm{norm}}^T \mathbf S^{\mathrm{raw}} \mathbf D_{\mathrm{norm}}
$$

对角元会变成 1，但非对角元一般仍然不是 0。

所以：

- **normalization** 只管“自己和自己”的范数
- **orthogonalization** 还要管“彼此之间”的重叠

这两个概念不能混在一起。

---

## 10. 和 p shell、d shell 的关系

归一化发生在径向收缩这一步，但一个 contracted radial function 最终还要乘上角向部分，才能变成 AO。

以 p shell 为例，一个 contracted radial function 可以写成：

$$
R_k(r_A)
=
\sum_{p=1}^{n_{\mathrm{prim}}} d_{pk}^{\mathrm{norm}} r_A^1 e^{-\alpha_p r_A^2}
$$

然后它会对应到三个角向分量：

$$
\boldsymbol\chi_{p,k}^{A}(\mathbf r)
=
R_k(r_A)
\begin{bmatrix}
\dfrac{x_A}{r_A} \\
\dfrac{y_A}{r_A} \\
\dfrac{z_A}{r_A}
\end{bmatrix}
$$

所以：

- `nctr = 1` 只表示“有 1 个 contracted radial function”
- 若是 p shell，它仍然会产生 3 个 p 型 AO
- 若是 spherical d shell，会产生 5 个 d 型 AO
- 若是 Cartesian d shell，会产生 6 个 d 型 AO

归一化不会改变 AO 个数，它改变的是“函数的尺度约定”。

---

## 11. 一个两 primitive 的最小例子

考虑一个 s shell（$l=0$），有两个 primitive：

$$
g_1(r)=e^{-\alpha_1 r^2},
\qquad
g_2(r)=e^{-\alpha_2 r^2}
$$

原始 contraction coefficient 为：

$$
\mathbf c
=
\begin{bmatrix}
c_1 \\
c_2
\end{bmatrix}
$$

### 第一步：primitive 归一化

$$
N_1 = \frac{1}{\sqrt{\mathrm{gaussian\_int}(2,2\alpha_1)}},
\qquad
N_2 = \frac{1}{\sqrt{\mathrm{gaussian\_int}(2,2\alpha_2)}}
$$

于是新的系数变成：

$$
\mathbf d
=
\begin{bmatrix}
N_1 c_1 \\
N_2 c_2
\end{bmatrix}
$$

### 第二步：构造 overlap 矩阵

$$
\mathbf S^{\mathrm{raw}}
=
\begin{bmatrix}
\mathrm{gaussian\_int}(2,2\alpha_1) & \mathrm{gaussian\_int}(2,\alpha_1+\alpha_2) \\
\mathrm{gaussian\_int}(2,\alpha_1+\alpha_2) & \mathrm{gaussian\_int}(2,2\alpha_2)
\end{bmatrix}
$$

### 第三步：contracted 归一化

$$
s = \frac{1}{\sqrt{\mathbf d^T \mathbf S^{\mathrm{raw}} \mathbf d}}
$$

最终存入 `libcint` 系数中的列向量是：

$$
\mathbf d^{\mathrm{norm}} = s\mathbf d
$$

也就是说：

$$
\chi(r) = d_1^{\mathrm{norm}} e^{-\alpha_1 r^2} + d_2^{\mathrm{norm}} e^{-\alpha_2 r^2}
$$

这里的 $d_1^{\mathrm{norm}}, d_2^{\mathrm{norm}}$ 已经同时吸收了：

- primitive normalization
- contraction normalization

---

## 12. `bas_ctr_coeff()` 返回的是什么

这是另一个很容易误解的点。

在 `_env` 中存下来的系数，本质上是：

$$
\mathbf C_{\mathrm{env}}
=
\mathbf N \, \mathbf C_{\mathrm{raw}} \, \mathbf S_c
$$

而 `bas_ctr_coeff()` 做了：

```python
cs = self._libcint_ctr_coeff(bas_id)
cs = numpy.einsum('pi,p->pi', cs, 1/gto_norm(l, es))
```

也就是说，它把 **primitive normalization** 从每一行除掉了：

$$
\mathbf C_{\mathrm{shown}}
=
\mathbf N^{-1} \mathbf C_{\mathrm{env}}
=
\mathbf C_{\mathrm{raw}} \, \mathbf S_c
$$

所以如果 `NORMALIZE_GTO=True`（PySCF 默认就是这样），那么 `bas_ctr_coeff()` 返回的并不是 basis 文件里的最原始系数，而是：

- 去掉了 primitive normalization
- 但仍然保留了 contracted normalization

因此你看到的数字可能和 basis 文件中的原始数字不完全一样。

---

## 13. 一句话总结

如果只记一条主线，可以记成：

```text
原始 basis coefficient
    -> 乘上 primitive normalization
    -> 再乘上每个 contraction 的整体归一化因子
    -> 存进 _env 给 libcint 直接使用
```

再压缩一点，就是：

```text
primitive normalization 负责“单个 primitive 的尺度”
contracted normalization 负责“线性组合后的整体范数”
```

而 PySCF 在 `make_bas_env()` 中选择的实现方式是：

```text
把这两层归一化都尽量吸收到 coefficient 矩阵里
```

这样后面的积分程序就可以直接读取 `_env` 中的系数，不需要再现场补这些归一化因子。
---





## 附录 A. 如何把 $\mathbf d_k^T \mathbf S^{\mathrm{raw}} \mathbf d_k$ 展开来看

如果你对这种矩阵写法不熟，最稳妥的读法是：**先把它翻译回双重求和，再翻译回积分。**

先记住三个对象：

$$
\mathbf d_k
=
\begin{bmatrix}
d_{1k} \\
d_{2k} \\
\vdots \\
d_{n_{\mathrm{prim}}k}
\end{bmatrix},
\qquad
\mathbf d_k^T
=
\begin{bmatrix}
d_{1k} & d_{2k} & \cdots & d_{n_{\mathrm{prim}}k}
\end{bmatrix}
$$

以及 overlap 矩阵：

$$
\mathbf S^{\mathrm{raw}}
=
\begin{bmatrix}
S_{11}^{\mathrm{raw}} & S_{12}^{\mathrm{raw}} & \cdots & S_{1n}^{\mathrm{raw}} \\
S_{21}^{\mathrm{raw}} & S_{22}^{\mathrm{raw}} & \cdots & S_{2n}^{\mathrm{raw}} \\
\vdots & \vdots & \ddots & \vdots \\
S_{n1}^{\mathrm{raw}} & S_{n2}^{\mathrm{raw}} & \cdots & S_{nn}^{\mathrm{raw}}
\end{bmatrix}
$$

这里为了简写，把 $n_{\mathrm{prim}}$ 暂时记成 $n$。

于是矩阵式

$$
\mathbf d_k^T \mathbf S^{\mathrm{raw}} \mathbf d_k
$$

完全展开就是：

$$
\sum_{p=1}^{n}\sum_{q=1}^{n}
d_{pk}\,S_{pq}^{\mathrm{raw}}\,d_{qk}
$$

也就是说：

- 左边那个 $d_{pk}$ 选中第 $p$ 个 primitive 的系数
- 中间的 $S_{pq}^{\mathrm{raw}}$ 表示第 $p$ 和第 $q$ 个 primitive 之间的 overlap
- 右边那个 $d_{qk}$ 选中第 $q$ 个 primitive 的系数
- 对所有 $p,q$ 都求和

把它再翻译回积分，就得到：

$$
R_k(r) = \sum_p d_{pk} g_p(r)
$$

因此

$$
\langle R_k | R_k \rangle
=
\left\langle \sum_p d_{pk} g_p \middle| \sum_q d_{qk} g_q \right\rangle
$$

利用积分的线性性，可以把求和提出来：

$$
\langle R_k | R_k \rangle
=
\sum_p \sum_q d_{pk} d_{qk} \langle g_p | g_q \rangle
$$

再定义

$$
S_{pq}^{\mathrm{raw}} = \langle g_p | g_q \rangle
$$

就得到

$$
\langle R_k | R_k \rangle
=
\sum_p \sum_q d_{pk} S_{pq}^{\mathrm{raw}} d_{qk}
=
\mathbf d_k^T \mathbf S^{\mathrm{raw}} \mathbf d_k
$$

所以这条公式本质上并不神秘，它只是把“双重求和”压缩写成了一个矩阵二次型。

最值得你记住的一个特例是两 primitive 情况。若

$$
\mathbf d_k =
\begin{bmatrix}
d_{1k} \\
d_{2k}
\end{bmatrix},
\qquad
\mathbf S^{\mathrm{raw}} =
\begin{bmatrix}
S_{11} & S_{12} \\
S_{21} & S_{22}
\end{bmatrix}
$$

那么

$$
\mathbf d_k^T \mathbf S^{\mathrm{raw}} \mathbf d_k
=
\begin{bmatrix}
d_{1k} & d_{2k}
\end{bmatrix}
\begin{bmatrix}
S_{11} & S_{12} \\
S_{21} & S_{22}
\end{bmatrix}
\begin{bmatrix}
d_{1k} \\
d_{2k}
\end{bmatrix}
$$

先算中间和右边：

$$
\mathbf S^{\mathrm{raw}} \mathbf d_k
=
\begin{bmatrix}
S_{11}d_{1k} + S_{12}d_{2k} \\
S_{21}d_{1k} + S_{22}d_{2k}
\end{bmatrix}
$$

再左乘 $\mathbf d_k^T$：

$$
\mathbf d_k^T \mathbf S^{\mathrm{raw}} \mathbf d_k
=
d_{1k}(S_{11}d_{1k} + S_{12}d_{2k}) + d_{2k}(S_{21}d_{1k} + S_{22}d_{2k})
$$

也就是

$$
\mathbf d_k^T \mathbf S^{\mathrm{raw}} \mathbf d_k
=
d_{1k}^2 S_{11} + d_{1k}d_{2k}S_{12} + d_{2k}d_{1k}S_{21} + d_{2k}^2 S_{22}
$$

如果 overlap 矩阵是对称的（通常确实如此，$S_{12}=S_{21}$），那么它会进一步化简成：

$$
\mathbf d_k^T \mathbf S^{\mathrm{raw}} \mathbf d_k
=
d_{1k}^2 S_{11} + 2d_{1k}d_{2k}S_{12} + d_{2k}^2 S_{22}
$$

这正对应你熟悉的“平方项 + 交叉项”结构。只不过在非正交基底下，交叉项的权重不是 1，而是 overlap $S_{12}$。

还有一个特别重要的类比：

- 如果 primitive 彼此正交归一，那么 $\mathbf S^{\mathrm{raw}} = \mathbf I$
- 这时
  $$
  \mathbf d_k^T \mathbf S^{\mathrm{raw}} \mathbf d_k = \mathbf d_k^T \mathbf d_k = \sum_p d_{pk}^2
  $$
- 也就是普通欧氏长度平方

所以你可以把

$$
\mathbf d_k^T \mathbf S^{\mathrm{raw}} \mathbf d_k
$$

理解成：**在非正交 primitive 基底里计算“长度平方”的正确方式。**
