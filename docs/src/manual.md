# 使用手册

本页按组件介绍核心 API 的用法与要点；逐项符号文档见 [API 参考](@ref)。目录：

1. [截断方案 `TruncationScheme`](@ref)
2. [格点 `ADTLattice` / `PTLattice`](@ref)
3. [系统算符 `AbstractImpurityOperator`](@ref)
4. [混合化样式 `HybridizationStyle`](@ref)
5. [影响泛函构建算法 `InfluenceFunctionalAlgorithm`](@ref)
6. [浴与关联函数（`ImpurityModelBase`）](@ref)
7. [ADT / PT 张量操作](@ref)
8. [观测量](@ref)
9. [MPO 哈密顿量（长程相互作用工具）](@ref)
10. [超参数与误差来源](@ref)
11. [代码结构](@ref)
12. [与文献的对应关系](@ref)

## 截断方案 `TruncationScheme`

张量网络计算中压缩键维的截断方案：

| 类型 | 构造方式 | 说明 |
|---|---|---|
| `TruncationDimCutoff` | `truncdimcutoff(D=χ, ϵ=ε, add_back=0)` | 同时限制最大键维 `D` 与截断阈值 `ϵ`（推荐） |
| `TruncateDim` | `truncdim(D)` | 只限制键维 |
| `TruncateCutoff` | `trunccutoff(ϵ=ε)` | 只按奇异值阈值截断 |
| `NoTruncation` | `NoTruncation()` | 不截断 |

预定义的默认值（`src/defaults.jl`）：

```julia
DefaultTruncation        # D=100, ϵ=1e-14
DefaultITruncation       # D=200,  ϵ=1e-10   （构建 IF 默认）
DefaultKTruncation       # D=1000, ϵ=1e-10   （系统动力学默认）
DefaultIntegrationTruncation  # D=10000, ϵ=1e-12
DefaultMPOTruncation     # D=10000, ϵ=1e-12
```

## 格点 `ADTLattice` / `PTLattice`

统一构造入口：

```julia
lattice = ADTLattice(N=Nt, δt=δt, contour=:real)      # 实时间
lattice = PTLattice(N=N,  δτ=δτ, contour=:imag)       # 虚时间
lattice = PTLattice(Nt=Nt, Nτ=Nτ, δt=δt, δτ=δτ, d=d, contour=:mixed)
```

常用属性与函数：

- `length(lattice)`：位点总数；`phydim(lattice)`：局域维数 `d`；
- `lattice.d`、`lattice.N` / `Nt` / `Nτ`、`lattice.δt` / `δτ`、`lattice.t`、`lattice.β`（虚时间时 `N*δτ`）；
- `branches(lattice)`：返回 `(:+, :-)`（实时间）或 `(:τ,)`（虚时间）等分支；
- `index(lattice, i, branch=:+)`：将 (时间步, 分支) 映射为位点位置；
- `ContourIndex(i, :+)`：轮廓索引对象；
- `vacuumstate(lattice)`：初始化真空态 ADT/PT。

!!! note "注意网格端点"
    `N` 个时间步覆盖的时间区间是 `t = N*δt`，采样点为 `i = 1..N` 对应 `t_i = i*δt`。与外部数据（如精确对角化）对比时注意网格是否包含端点、起始点，建议按时间显式对齐而不是按下标对齐。

## 系统算符 `AbstractImpurityOperator`

| 类型 | 构造 | 说明 |
|---|---|---|
| `ImpurityHamiltonian` | `ImpurityHamiltonian(m::Matrix)` | 幺正演化 $e^{\mp i \hat{H}_S \delta t}$ |
| `ImpurityLindbladian` | `ImpurityLindbladian(H, jumpops)` | 含耗散（Lindblad）演化，`jumpops` 为跃迁算符列表 |

由 `sysdynamics(lattice, model, trunc=...)` 生成裸系统动力学张量 `mpsK`（ADT 或 PT）。

## 混合化样式 `HybridizationStyle`

描述系统-浴耦合的具体形式：

| 类型 | 耦合形式 | 约束 | 框架 |
|---|---|---|---|
| `AdditiveHyb(op)` | $\hat{A}(\hat{b}^\dagger+\hat{b})$ 型，对角 | `op` 为向量或对角矩阵（厄米） | ADT |
| `NonAdditiveHyb(op)` | $\hat{A}(\hat{a}+\hat{a}^\dagger)$，$\hat{A}=\hat{A}^\dagger$ | `op` 厄米矩阵 | PT |
| `NonDiagonalHyb(op)` | $\hat{A}\hat{a} + \hat{A}^\dagger \hat{a}^\dagger$ | `op` 任意方阵 | PT（文献的核心情形） |

含时版本（`TdHybridizationStyle`）：`AdditiveTdHyb`、`NonAdditiveTdHyb`、`NonDiagonalTdHyb`。

## 影响泛函构建算法 `InfluenceFunctionalAlgorithm`

```julia
# 1) 部分影响泛函（Partial IF）：仅适用于对角耦合 AdditiveHyb
alg = PartialIF(trunc=trunc)

# 2) 平移不变影响泛函（XTRG-IF，XTRG 式），文献推荐
alg = XTRGIF(;
    algexpan = OverDeterminedProny(n=20, tol=1.0e-8, verbosity=2),  # 混合化函数指数展开
    algevo   = WII(),        # 或 WI()、ComplexStepper()、FirstOrderStepper()
    algmult  = DMRGMult1(trunc, initguess=:rand, maxiter=10),   # 或 SVDCompression(trunc)
    k        = 7,            # XTRG 步数：时间步长 1/2^k
    fast     = true,         # true：树形二分方案（约 k 次乘法）；false：串行 2^k-1 次
    verbosity= 0,
)
```

- `algexpan`：`ExponentialExpansionAlgorithm`，包括 `OverDeterminedProny`、`DeterminedProny`（`exponential_expansion`、`expansion_error` 可用于误差分析，由重导出的 `ExpExp` 包提供）；
- `algevo`：`TimeEvoMPOAlgorithm`，把 $\hat{H}_{\text{eff}}$ 指数化为 MPO 的步进器（`WI`/`WII`/`FirstOrderStepper`/`ComplexStepper`）；
- `algmult`：`DMRGAlgorithm`，MPO-MPO 乘法的压缩算法，`DMRGMult1`（单站点 DMRG 迭代，`initguess ∈ {:svd, :pre, :rand}`，`maxiter`）或 `SVDCompression`。

构建 IF 的入口函数：

```julia
mpsI = hybriddynamics(lattice, corr, hyb)                 # 默认 PartialIF / 默认截断
mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)    # 对角耦合 + 部分 IF
mpsI = hybriddynamics(lattice, corr, hyb, alg)            # 指定 IF 算法（对角或非对角）
mpsI = hybriddynamics!(gmps, lattice, corr, hyb, alg)     # 就地版本
mpsI = hybriddynamics_naive(lattice, corr, hyb, trunc=trunc)  # N² 门操作的朴素参考实现（对角）
```

> 对角耦合，`hybriddynamics` 既可用 `PartialIF`（逐位点乘，每因子键维为 2），也可用 `XTRGIF`。对非对角耦合只能使用 `XTRGIF`。

## 浴与关联函数（`ImpurityModelBase`）

```julia
spec = spectrum(f, lb=0, ub=wc)          # 由函数定义谱密度 J(ω)
spec = Leggett(d=1, ωc=5, α=0.1)         # 预定义 Leggett 谱
bath = bosonicbath(spec, β=β)            # 玻色浴（有限温度 β）
corr = correlationfunction(bath, lattice) # 离散化到格点上的关联函数 Δ(τ,τ')
```

`correlationfunction` 会根据格点轮廓自动调用 `Δt`（实时间）、`Δτ`（虚时间）或 `Δm`（混合轮廓）。

要点：
- 零温取 `β=Inf`；
- `spectrum` 的积分区间 `[lb, ub]` 应覆盖谱的主要权重（如 Ohmic 谱取 `ub=3~5ωc`）；若 $J(\omega)$ 在 $\omega=0$ 发散或与发散核卷积，需从非零下界开始（详见[实践指南](@ref)）。

## ADT / PT 张量操作

- `mult(a, b, algmult)` / `mult!(a, b, trunc)`：张量网络乘法（元素级乘积 / MPO 乘法）；
- `integrate(lattice, args...)` / `integrate(mpsA, mpsB)`：求配分函数（路径积分求和）；
- `apply!(term, mps)`：应用局域算符；
- `canonicalize!`、`leftorth!`、`rightorth!`：正交化（可指定 `Orthogonalize`）；
- `bond_dimension(mps)`、`bond_dimensions(mps)`：键维查询；
- `distance(mps1, mps2)` / `distance2`：两个张量的距离（相对误差验证用）；
- `randomadt` / `randompt`：随机张量（测试用）；`vacuumstate(lattice)`：真空态。

## 观测量

共有两条测量路径，均适用于 ADT 与 PT：

**路径 A：算符插入（任意算符，含非对角与两点关联）**

把算符（可以是任意矩阵）作为 `ContourOperator` 插入系统动力学，再整体求值：

```julia
# 单点：⟨op(t_i)⟩
ct   = ContourOperator(ContourIndex(i), op)
mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)    # 实时间需 ρ₀；虚时间不用
mps2 = mult!(mpsK, mpsI, trunc=trunc)
v    = integrate(mps2) / integrate(mps)

# 两点关联：⟨op2(t_i) op1(t_j)⟩
ct   = ContourOperator([ContourIndex(i), ContourIndex(j)], [op2, op1])
mpsK = sysdynamics(lattice, model, ct, trunc=trunc)
# ……同样 boundarycondition! → mult! → integrate/integrate
```

**路径 B：环境缓存（批量单点观测量）**

先建一次缓存，再逐点测（已归一化）：

```julia
# PT 框架（非对角耦合，或需要非对角观测量）
cache = environments(lattice, mps, ρ₀=ρimp)
v = expectationvalue(ContourOperator(ContourIndex(i, :+), op), cache)

# ADT 框架（对角耦合，对角观测量）
cache = environments(lattice, mpsK, mpsI)
v = expectationvalue(ADTTerm(index(lattice, i, branch=:+), zdiag), cache)
```

虚时间 / 混合轮廓同样支持（PT 用 `environments(lattice, mps)`，ADT 用 `environments(lattice, mpsK, mpsI)`）。ADT 的多点形式 `ADTTerm((pos2, pos1), (v2, v1))` 可测对角两点关联。

辅助函数：`Zvalue(cache)`（配分函数）、`Zvalue2(cache)`、`TransferMatrix`（转移矩阵）、
`correlation(lattice, model, op, mpsI[, ρ0])`（两点关联函数）、`heatcurrents`。

## MPO 哈密顿量（长程相互作用工具）

`src/mpohamiltonian/` 提供了构造长程（指数衰减）相互作用哈密顿量 MPO 的独立工具，用于其他一维量子多体问题：

```julia
# SchurMPOTensor：把 [局域项 + 指数衰减长程项] 编码为紧凑 MPO 站点张量
h = SchurMPOTensor(h1, h2s)    # h2s 为 ExponentialDecayTerm / GenericDecayTerm / PowerlawDecayTerm 列表
mpo = MPOHamiltonian([h for _ in 1:L])
tensors = tompotensors(mpo)              # 转成稠密 MPO 站点张量
tensors2 = timeevompo(tensors, dt, WII())   # 时间演化（WI / WII / ComplexStepper / FirstOrderStepper）
```

## 超参数与误差来源

本方法共有四类可控误差来源（文献 Sec. IV）：

| 超参数 | 含义 | 建议默认 | 误差来源 |
|---|---|---|---|
| `δt` / `δτ` | 轮廓离散步长 | 收敛性检查确定 | QuAPI 一阶 Trotter 分解误差，整体 $O(t\delta t)$ |
| `χ`（`trunc.D`） | MPO/MPS 键维截断 | 文献示例 30–200 | MPO 压缩（截断）误差 |
| `n`（`algexpan`） | 混合化函数指数展开项数 | 文献默认 20 | Prony 近似误差 |
| `m`（`alg.k`） | XTRG 步数（时间步 $1/2^m$） | 文献默认 7 | XTRG 热态构造误差（指数收敛） |
| `d`（格点） | 玻色杂质局域 Hilbert 空间截断 | 依温度而定 | 局域截断误差 |

要点（均已在文献中以 JC 模型、双单模浴、非相互作用玻色模、sub-Ohmic 浴等算例验证）：

- 误差随 `χ`、`m`、`n` 增加快速饱和；`δt` 越小误差单调下降；
- 强耦合一般比弱耦合更难收敛（需要更大 `χ`）；
- JC 型（粒子数守恒）耦合生成的纠缠远少于 Rabi 型耦合，通常收敛更快；
- 零温情形（$\beta=\infty$）一般比有限温更容易收敛；
- 计算成本主要来自 XTRG 中的 MPO-MPO 乘法，理论标度为 $O(N\chi^4 d^3)$。

## 代码结构

```
src/
├── TEMPO.jl                    # 模块定义与导出符号
├── tensorops/                  # 截断、张量操作、正交化等基础工具（algorithms.jl 为 DMRG 乘法/MPS 算法）
├── defaults.jl                 # 默认超参数
├── mpohamiltonian/             # MPO 哈密顿量（SchurMPO、长程项、时间演化步进器）
├── contourindices.jl           # ContourIndex、branch
├── adt/  pt/                   # 增强密度张量（MPS）与过程张量（MPO）数据结构及运算
├── conversions.jl              # PT ↔ ADT 转换辅助（copy 张量等）
├── adtterms.jl / fockterms.jl  # ADTTerm / FockTerm / ProdFockTerm 等局域项
├── adtlattices/ ptlattices/    # 实/虚/混合轮廓格点定义（Fock 排序）
├── contouroperators.jl         # ContourOperator（PT 上的算符）
├── correlationfunction.jl      # 浴关联函数到格点的离散化
├── influencefunctional/        # Feynman-Vernon IF：PartialIF / 平移不变 IF（ADT 与 PT 两种）
├── tdinfluencefunctional/      # 含时耦合的 IF 与含时混合化样式
├── boundarycondition.jl        # 初始态 / 边界条件
├── models/                     # 幺正（ImpurityHamiltonian）与耗散（ImpurityLindbladian）动力学
└── observables/                # 环境缓存、期望值、转移矩阵、关联函数、热流

docs/tutorials/                  # 教程：spinboson 与三篇论文复现 notebook
benchmark/                       # 基准算例：自旋玻色子、独立玻色模、实时间等
test/                            # 测试套件（含与 Exact Diagonalization 的对照）
```

## 与文献的对应关系

| 文献内容 | 工具包实现 |
|---|---|
| 对角 + 对易耦合（标准 TEMPO） | `ADTLattice` + `AdditiveHyb` + `PartialIF`/`XTRGIF` |
| 对角 + 非对易耦合（多浴） | 多个 `AdditiveHyb` 的 IF 相乘（`mult!`） |
| 非对角耦合（共轭对，核心推广） | `PTLattice` + `NonDiagonalHyb` + `XTRGIF` |
| 实时间 Keldysh 轮廓 | `contour=:real` |
| 虚时间轮廓 | `contour=:imag` |
| L 形 Kadanoff–Baym 轮廓 | `contour=:mixed`（`MixedPTLattice`/`MixedADTLattice`） |
| QuAPI 离散化（附录 C） | `correlationfunction(bath, lattice)` 中的 `Δt`/`Δτ`/`Δm` |
| 混合化函数的指数展开（附录 D） | `OverDeterminedProny`/`DeterminedProny`（`algexpan`） |
| XTRG 构造有效热态（文献 Fig. 3） | `XTRGIF(k=..., fast=...)` |
| PT 中系统哈密顿量的吸收（文献 Fig. 2b） | `sysdynamics(lattice, model)` + `mult!` |
| 观测量计算（文献 Fig. 1d, 1e） | `environments` + `expectationvalue` |
| JC 自旋玻色子模型 | `docs/tutorials/spinboson/jctype.jl` |
| 标准自旋玻色子模型 | `docs/tutorials/spinboson/rabitype.jl`、`benchmark/sb.jl` |
| 非相互作用/相互作用玻色杂质 | `benchmark/independentbosons.jl`、`bosonicimpurity.jl` |
