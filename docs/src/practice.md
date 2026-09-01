# 实践指南

本页总结在复现文献算例（`docs/tutorials/` 下的 spinboson、strathearn2018、otterpohl2025、guo2026）过程中积累的实战经验：如何选择计算路径、如何做收敛检验、常见的坑以及如何诊断错误结果。**建议在做正式计算之前通读一遍。**

## 1. 选择计算路径

先回答两个问题：

1. **耦合算符是什么形式？** 对角（厄米、与 $H_S$ 对易）→ ADT 框架；非对角（如 $\sigma_-$）或非对易多浴 → PT 框架。
2. **要测什么观测量？** 只测对角量（如 $\langle S_z(t)\rangle$）→ ADT 够用；要测非对角量（如 $\langle\sigma_x(t)\rangle$、$\langle a\rangle$）→ 必须走 PT 路径。

| 情形 | 推荐路径 | 理由 |
|---|---|---|
| 对角耦合 + 对角观测量 + d=2 | `ADTLattice` + `AdditiveHyb` + `PartialIF`（默认） | 最快，每因子键维 2 |
| 对角耦合 + 对角观测量 + d≥3（如 spin-1） | `ADTLattice` + `AdditiveHyb` + `TranslationInvariantIF` | PartialIF 因子数随 d 增长，TTI-IF 明显更快 |
| 需要非对角观测量 | `PTLattice` + `NonDiagonalHyb` + `TranslationInvariantIF`；或 ADT 路径 + 算符插入（见 §3） | PT 路径批量测更高效；ADT 算符插入适合少量算符/两点关联 |
| 非对角/非厄米耦合 | `PTLattice` + `NonDiagonalHyb` + `TranslationInvariantIF` | 唯一支持 |

!!! tip "对角耦合也能测非对角观测量"
    两条路线：① 对**实**对角算符 `op`，`NonDiagonalHyb(op)` 与 `AdditiveHyb(op)` 物理等价，但走 PT 路径、可用 `ContourOperator` 批量测非对角观测量（如 otterpohl2025 教程中 `NonDiagonalHyb(Matrix(Diagonal(0.5 .* [1.0, -1.0])))` 即 $\sigma_z/2$ 耦合，同时测 $P=\langle\sigma_z\rangle$ 与 $C=\langle\sigma_x\rangle$）；② 保持在 ADT 路径，用**算符插入**（`sysdynamics(lattice, model, ct)`，见 §3），适合少量算符或两点关联函数。

## 2. 收敛检验：流程与经验值

四类误差来自 `δt`、`χ`、`n`（指数展开项数）、`k`（XTRG 步数）。推荐流程：

1. 固定较大 `χ`、`n=20`、`k=7`，扫描 `δt`（如 0.2 → 0.1 → 0.05），取结果不再变化的 `δt`；
2. 固定该 `δt`，扫描 `χ`（如 25 → 50 → 100），直到观测量（尤其是峰值、极小值位置）不再变化；
3. 若用 TTI-IF，检查 `exponential_expansion` 的展开误差（见 §5）。

实测经验值（可作起点，仍需按自己的参数验证）：

| 算例 | 难点 | 经验值 |
|---|---|---|
| Ohmic SBM（$\Omega=1$，$\omega_c=5$，$T=0$）**穿过局域化相变** | 相变点附近纠缠极大 | $\chi\ge 50$ 且 $\delta t\le 0.1$（$\chi=25$、$\delta t=0.15$ 时 $S_z$ 误差可达 ~0.04） |
| Ohmic SBM 低温 | — | $\beta=20$ 已等效零温（与 $\beta=100$ 差别 <0.001） |
| $1/f^\eta$ 噪声（$s<0$，Otterpohl 型） | 时间幂律记忆使 IF 键维随时间持续增长 | 短时间（$t\le5$）$\chi=40$ 够用；较强耦合 $\alpha=0.04$ 需 $\chi=100$ |
| JC 型（粒子数守恒）非对角耦合 | 纠缠远少于 Rabi 型 | $\chi=30$ 通常即可 |
| spin-1（d=3）+ TTI-IF | 局域维数大 | $\chi=50$、$\delta t=0.1$、$n=20$、$k=7$ |

## 3. 观测量构造：最容易踩的坑

### 混合化样式的算符参数

- `AdditiveHyb` 接受**对角元向量**或对角矩阵。传矩阵 `z` 时会被取对角元，但为清晰起见建议直接传 `zdiag = diag(z)`；
- `NonDiagonalHyb(op)` 的含义是 $\hat{A}\hat{a} + \hat{A}^\dagger\hat{a}^\dagger$，即 `op` 是耦合算符 $\hat{A}$ 本身（**含系数**）。自旋玻色子 $\tfrac{1}{2}\sigma_z\hat\xi$ 应写 `NonDiagonalHyb(Matrix(Diagonal(0.5 .* [1.0, -1.0])))`，即 $\hat{A}=\sigma_z/2$；
- 注意文献间的算符约定差异：有的论文用 $S=\sigma/2$（TEMPO 包约定，Pauli 算符无两倍因子），对比论文数据时先核对耦合算符与 $H_S$ 的定义（如 $H_s=\Omega\sigma_x/2$ 对应 `ImpurityHamiltonian(Ω .* Sx/2)` 还是 `Ω .* Sx`，取决于 $S$ 还是 $\sigma$）。

### ADT 路径：`ADTTerm` 是单点对角算符，非对角量用算符插入

```julia
m = ADTTerm(index(lattice, i, branch=:+), zdiag)   # zdiag 是对角元向量
v = expectationvalue(m, cache)                     # 已归一化
```

- `ADTTerm(pos, data)` 中 `data` 是**对角元向量**，作用在单个位点上（系统 + 分支配对位点已被格点吸收为一步）；
- 多位置形式 `ADTTerm((pos2, pos1), (v2, v1))` 可测**对角**两点关联：`mps2 = apply!(m, copy(mps)); v = integrate(mps2) / integrate(mps)`；
- **不要**把完整矩阵传给单点 `ADTTerm`；
- **非对角单点算符或任意两点关联**用算符插入（对 ADT 与 PT 都适用，见 `test/models/rabimodel.jl`）：

  ```julia
  ct   = ContourOperator([ContourIndex(i), ContourIndex(1)], [op2, op1])  # op 可为任意矩阵
  mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
  mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)   # 实时间需 ρ₀，虚时间不用
  mps2 = mult!(mpsK, mpsI, trunc=trunc)
  v    = integrate(mps2) / integrate(mps)             # 已归一化
  ```

  每次插入都要重建 `mpsK` 并重乘 IF，适合少量算符；批量单点测量优先用环境缓存路径。

### PT 路径：`ContourOperator` 接受任意矩阵

```julia
cache = environments(lattice, mps, ρ₀=ρimp)        # 实时间 PT 必须给 ρ₀！
m = ContourOperator(ContourIndex(i, :+), x)        # x 是任意矩阵
v = expectationvalue(m, cache)
```

- 实时间 PT 的 `environments` 需要 `ρ₀`（初态密度矩阵），漏掉会得到默认最大混合态的结果；
- 两点插入用 `ContourOperator([c1, c2], [op1, op2])` 并通过 `sysdynamics(lattice, model, ct)` 进入系统动力学（虚时间 Green 函数的标准做法）。

### 初态

```julia
mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)  # ADT 路径
```

- `ρ₀` 可以是密度矩阵，也可以是对角向量（纯态对角元）；
- 初态施加在 `mpsK`（系统动力学）上，不要施加在 `mpsI`（影响泛函）上。

## 4. 结果非物理？先怀疑 χ 不足

χ 不足的**典型症状不是报错，而是静默的错误结果**：

- 归一化漂移：$P(t)>1$、$|C(t)|>1$；
- 出现非物理振荡或极小值位置明显偏移；
- 长时间行为发散。

例如 otterpohl2025 教程中 $\alpha=0.04$、$\chi=40$ 时得到 $P(\text{end})=2.94$、$C_{\min}=-1.58$（非物理），换 $\chi=100$ 后恢复 $P\le1$、$C$ 出现物理振荡。**任何新算例都应至少做两个 χ 值的对比。**

## 5. 谱密度与关联函数

- `spectrum(f, lb=0, ub=...)`：上限应取到谱权重基本衰减完（Ohmic/sub-Ohmic 取 `ub=3~5ωc` 足够）；
- **低频发散谱**：若 $J(0)\neq 0$ 且与玻色因子 $1/\epsilon$ 卷积，$\omega=0$ 处对数发散。处理办法是从小的非零下界积分（如 segal2023 教程从 $\gamma/100$ 开始）；
- **零温**：`bosonicbath(spec, β=Inf)`；
- **`ExpExp` 警告** `can not find a good approximation with L=..., n=...`：表示 Prony 展开未达到 `tol`。对 $s<0$ 的幂律记忆核这很常见，通常仍可计算（TEMPO 只需关联函数），但应：
  - 增大 `n`（如 20 → 30）或放宽 `tol`（如 `1e-8` → `1e-6`）并对比结果；
  - 警告只在收敛性检查时认真对待，最终结果应用无警告（或误差可忽略）的参数复核。

## 6. 性能：缓存与算法选择

- **缓存影响泛函 MPS/MPO**：TTI-IF 的 XTRG 构造可能耗时小时级（如 $\chi=100$、$N=103$ 约 3.7 小时）。建议用 `Serialization` 把 `mpsI` 存盘，重跑观测量时直接反序列化：

  ```julia
  using Serialization
  mpspath = "data/mypi_beta20_dt0.1_alpha0.3_chi50_N200.mps"
  mpsI = ispath(mpspath) ? Serialization.deserialize(mpspath) :
         (I = hybriddynamics(lattice, corr, hyb, alg); Serialization.serialize(mpspath, I); I)
  ```
- **TTI-IF 开两个开关**：`fast=true`（树形二分，约 `k` 次乘法而非 $2^k-1$ 次）、`k` 适当（`k=7` 即最小步 $1/128$，更大 `k` 收益递减）；
- **PartialIF 对 d=2 是最快的**，不要为了"新算法"盲目换 TTI-IF；d≥3 时 TTI-IF 优势明显；
- 设置 `OMP_NUM_THREADS=1`（或按机器适当开线程）避免 BLAS 过订阅；DMRG 乘法本身是单线程的。

## 7. 与参考数据（论文/ED）对比

1. **先核对约定**：耦合算符、$H_S$、谱密度的归一化（如 $J(\omega)=2\alpha\omega e^{-\omega/\omega_c}$ 是否带 $\pi$）、时间单位；
2. **按时间对齐网格**：TEMPO 的采样点是 `t_i = i*δt, i=1..N`，与外部数据（ED 网格、数字化数据）对比时用时间值索引而不是下标：
   ```julia
   n_on_t = [ne[round(Int, t_i / 0.005) + 1] for t_i in t]   # ED 步长 0.005
   ```
3. **数字化论文数据**：注意图内 legend 与论文文字描述的一致性（strathearn2018 教程中曾因误用另一 panel 的 legend 把 α 标错）；数字化曲线自带噪声，对比关注趋势与极值位置而不是逐点差；
4. **物理量纲**：浴紫外截断 $1/\omega_c$ 决定赝相干区极小值位置等特征时间尺度，量纲要自洽。

## 8. 快速诊断表

| 症状 | 可能原因 | 处理 |
|---|---|---|
| $P(t)>1$、$\lvert C\rvert>1$ | χ 不足 | 增大 χ 重跑（§4） |
| 结果随 t 发散 | χ 不足或 δt 过大 | 先加倍 χ，再减半 δt |
| 极小值/峰值位置系统性偏移 | δt 离散误差 | 扫描 δt（§2） |
| `ExpExp` 不收敛警告 | 关联函数难展开 | 增大 n 或放宽 tol，复核（§5） |
| 期望值对应最大混合态 | PT 路径忘传 `ρ₀` | `environments(lattice, mps, ρ₀=ρimp)`（§3） |
| `expectationvalue` 类型错误 | ADT/PT 项混用 | `ADTTerm` 配 ADT 缓存，`ContourOperator` 配 PT 缓存（§3） |
| 浴关联函数发散 | 谱低频发散 | 谱积分下界取小正值（§5） |
| 与论文对不上 | 算符/谱约定差异 | 核对 $S$ vs $\sigma$、$\pi$ 归一、单位（§7） |
