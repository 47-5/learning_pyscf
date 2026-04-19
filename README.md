# learning_pyscf

这个仓库只做三件事：

1. 按周推进 PySCF 源码阅读和理论学习。
2. 在 `annotated_pyscf_source/` 里保存我正在读或已经读过的 PySCF 源码副本，并加入自己的注释。
3. 在 `notebooks/` 里逐步放少量可复现实验 notebook。

## 使用原则

- 不追求一开始就很完整，先把 RHF/RKS 主线读透。
- 注释源码时保留原文件相对路径，例如 `pyscf/scf/hf.py` 可以复制到 `annotated_pyscf_source/pyscf/scf/hf.py`。
- 每次只注释一个小片段：一个函数、一个类、或一次 SCF/DFT 数据流。
- 不上传 Gaussian、ORCA、CP2K 的大输出文件或未发表科研体系。

## 24 周详细计划表

| 周 | 主题 | 理论阅读 | PySCF 源码阅读 | 对照经验 | 本周产出 |
|---|---|---|---|---|---|
| 1 | 仓库和学习地图 | 复习 PySCF User Guide 的分子输入；明确 charge、spin、basis、unit 的含义 | 先看 `examples/gto/00-input_mole.py`，浏览 `pyscf/gto/mole.py` 的 `Mole` 对象 | 对照 Gaussian/ORCA 的坐标、电荷、多重度输入 | 在 README 中勾画个人目标；复制并注释第一个小例子 |
| 2 | 分子对象和基组 | 读 Szabo & Ostlund 中 AO/MO、基函数、Slater determinant 的基础部分 | 读 `examples/gto/04-input_basis.py`、`examples/gto/11-basis_info.py` | 对照 Gaussian/ORCA 中 STO-3G、def2-SVP、cc-pVDZ 的写法 | 注释一个基组输入例子；记录 H2O 不同基组的 AO 数 |
| 3 | 一电子积分 | 复习 overlap、kinetic、nuclear attraction integrals | 读 `examples/gto/20-ao_integrals.py`，再看 `pyscf/scf/hf.py` 的 `get_hcore`、`get_ovlp` | 对照 Gaussian/ORCA 输出里 basis functions、one-electron energy 的概念 | 注释 `get_hcore` 和 `get_ovlp`；写清楚 `S`、`T`、`V` 的矩阵维度 |
| 4 | 二电子积分和 ERI | 复习 `(mu nu|lambda sigma)`、Coulomb/exchange 的来源 | 继续读 `examples/gto/20-ao_integrals.py`，只浏览 `pyscf/scf/_vhf.py` 入口 | 理解为什么大体系中 HF/hybrid 会变慢 | 注释 ERI 示例；画出 ERI 四指标张量到 J/K 的关系 |
| 5 | RHF 主循环总览 | 读 Roothaan-Hall 方程 `FC = SCE` | 读 `pyscf/scf/hf.py` 的 `kernel` 主循环 | 对照 Gaussian 的 `SCF Done` 和 ORCA 的 SCF iteration 表 | 注释 `kernel` 的循环结构；画出 `D -> J/K -> F -> C -> D` |
| 6 | 密度矩阵 | 复习 closed-shell density matrix 的构造 | 读 `make_rdm1`、`get_occ`、`eig` | 对照轨道占据数、HOMO/LUMO 输出 | 注释 `make_rdm1`；解释 `mo_coeff`、`mo_occ`、`dm` |
| 7 | Fock 矩阵和能量 | 读 HF 电子能表达式，区分 `E1`、`Ecoul`、`Etotal` | 读 `get_fock`、`energy_elec` | 对照 Gaussian/ORCA 总能量分量 | 注释 `energy_elec`；写一页公式到变量名对照 |
| 8 | J/K 构造 | 复习 Coulomb matrix 和 exchange matrix | 读 `get_jk`、`get_veff`，只跟到 Python 层接口 | 对照 hybrid functional 为什么比 GGA 慢 | 注释 `get_veff`；说明 `vj`、`vk`、`vhf` 的关系 |
| 9 | SCF 收敛和 DIIS | 读 DIIS 的误差向量思想，不追求推导全部细节 | 读 `pyscf/scf/diis.py` 和 `pyscf/lib/diis.py` 的主接口 | 对照 Gaussian/ORCA 的 DIIS、MaxCycle、level shift | 注释 DIIS error vector；记录一个不收敛时可能调什么 |
| 10 | UHF/开壳层入门 | 复习 alpha/beta 密度矩阵和 spin contamination 概念 | 浏览 `pyscf/scf/uhf.py`，重点看和 RHF 不同的数据形状 | 对照 ORCA/Gaussian 的 multiplicity 和 unrestricted 关键词 | 注释一个 UHF 小片段；写 RHF/UHF 数据结构差异表 |
| 11 | RKS 总览 | 读 Kohn-Sham DFT 基本方程；理解和 HF 共用 SCF 框架 | 读 `pyscf/dft/rks.py` 的 `RKS`、`KohnShamDFT`、`get_veff` | 对照 Gaussian/ORCA 的 PBE、B3LYP、PBE0 关键词 | 注释 RKS 类结构；说明为什么 DFT 仍然是 SCF |
| 12 | 数值积分网格 | 读 Becke partition、radial grid、angular grid 的基本概念 | 读 `pyscf/dft/gen_grid.py` 的 `Grids` 和 `gen_atomic_grids` | 对照 Gaussian `Integral=UltraFine`、ORCA `GridX` | 注释 `Grids`；记录 grid level 对应的是数值精度问题 |
| 13 | AO 到 rho(r) | 读 `rho(r) = sum D_mn chi_m(r) chi_n(r)` | 读 `pyscf/dft/numint.py` 的 `eval_ao`、`eval_rho` | 对照 CP2K 中 density grid 的直觉 | 注释 `eval_ao`、`eval_rho`；画出 `dm + ao -> rho` |
| 14 | XC functional 调用 | 读 LDA/GGA/meta-GGA 需要的变量：rho、grad rho、tau | 读 `pyscf/dft/numint.py` 的 `eval_xc` 相关入口，浏览 `pyscf/dft/libxc.py` | 对照 Libxc、Gaussian/ORCA 泛函命名差异 | 注释一次 PBE 或 LDA 的 XC 调用路径 |
| 15 | RKS 数值积分回 AO 矩阵 | 读 `v_xc(r)` 如何投影回 AO basis 的想法 | 读 `numint.nr_rks`，重点看 block loop 和 `vmat` 的构造 | 对照为什么 grid 越密越慢 | 注释 `nr_rks` 的大块结构；写 `rho -> exc/vxc -> vmat` 流程 |
| 16 | Hybrid DFT | 复习 exact exchange、hybrid functional 的组成 | 回看 `rks.get_veff` 中 hybrid 分支和 `get_jk` 调用 | 对照 B3LYP/PBE0 比 PBE 慢的实际经验 | 注释 hybrid 分支；总结 GGA 和 hybrid 在代码路径上的差异 |
| 17 | DFT 网格收敛实验 | 读数值积分误差和能量收敛的关系 | 不新增大段源码，回看 `Grids` 和 `nr_rks` | 对照 Gaussian/ORCA 不同 grid 设置结果可能不同 | 在 `notebooks/` 中建立第一个 notebook：H2O/PBE grid level 扫描 |
| 18 | Density fitting/RI 概念 | 读四中心积分到三中心积分的 RI 近似思想 | 浏览 `pyscf/df/df.py`、`pyscf/df/df_jk.py` 的入口 | 对照 ORCA 的 RI、RIJCOSX、AutoAux、def2/J | 注释 DF 入口；写普通 SCF vs density fitting 的差异 |
| 19 | SCF callback 和可观测量 | 读 SCF 每轮应记录哪些量：E、dE、ddm、gradient norm | 回看 `kernel` 中 callback 触发位置 | 对照 Gaussian/ORCA 输出表格 | 建一个 notebook 或小脚本，记录每轮 SCF 收敛轨迹 |
| 20 | Broken symmetry 和稳定性 | 读稳定性分析、开壳层初猜的概念 | 浏览 `examples/scf/17-stability.py` 和 `examples/dft/32-broken_symmetry_dft.py` | 对照过渡金属/自由基计算中的 broken symmetry 经验 | 注释一个 broken-symmetry 示例；记录适用场景和风险 |
| 21 | PBC Cell 入门 | 复习晶胞、周期边界、赝势、k 点的概念 | 读 `examples/pbc/00-input_cell.py` 和 `pyscf/pbc/gto/cell.py` 入口 | 对照 CP2K 的 `CELL`、`COORD`、`BASIS_SET`、`POTENTIAL` | 注释一个 `Cell` 输入例子；建立 CP2K 到 PySCF PBC 对照表 |
| 22 | PBC Gamma-point SCF/DFT | 读 gamma point 近似和分子计算的相似处 | 读 `examples/pbc/10-gamma_point_scf.py`、浏览 `pyscf/pbc/scf/hf.py` | 对照 CP2K 单胞 gamma point 计算 | 注释 gamma-point SCF 示例；说明 `Mole` 和 `Cell` 的差异 |
| 23 | k 点采样 | 读 Bloch theorem 和 k-point sampling 的实用意义 | 读 `examples/pbc/20-k_points_scf.py`，浏览 `pyscf/pbc/scf/khf.py` | 对照 CP2K `KPOINTS` | 注释 k-point 示例；记录 gamma 与 k-point 数据结构差异 |
| 24 | 阶段总结 | 回看 HF、RKS、grid、DF、PBC 之间的主线 | 回顾自己注释过的源码文件 | 选择下一阶段方向：XC 实现、PBC、性能、响应性质或贡献文档 | 写一份 2-3 页总结：我如何从软件使用者走向源码理解 |

## 每周最小执行模板

| 时间 | 任务 |
|---|---|
| 第 1 次 1 小时 | 读理论 20 分钟；读源码 40 分钟 |
| 第 2 次 1 小时 | 继续读源码并加注释 |
| 第 3 次 1 小时 | 跑一个最小例子或 notebook |
| 第 4 次 1 小时 | 整理本周理解、疑问和下一周入口 |

## 注释源码的建议格式

```python
# MYNOTE: 这里对应 Roothaan-Hall 方程中的 F 矩阵构造。
# MYNOTE: dm 是 AO basis 下的 density matrix，shape 通常是 (nao, nao)。
```

建议只用 `MYNOTE:` 标记自己的注释，这样以后可以全文搜索，也方便和原始 PySCF 源码区分。