# TEMPO.jl 工具包文档

本工具包是 **时间演化矩阵乘积算符（Time-Evolving Matrix Product Operator, TEMPO）** 方法的一个 Julia 实现，其算法理论基于文献：

> C. Guo, W. Wu, X. Xu, T. Jiang, P.-X. Chen, and R. Chen,
> *Time-evolving matrix product operators for off-diagonal system-bath coupling*,
> **Phys. Rev. B 114, 125413 (2026)**（即 `docs/GuoChen2026a.pdf`）。

与只支持对角系统-浴耦合的原始 TEMPO 不同，本实现基于**过程张量（Process Tensor, PT）**框架，将 TEMPO 推广到了更一般的**非对角系统-浴耦合**情形（系统通过一对共轭的非厄米算符 $\hat{A}^\dagger, \hat{A}$ 与浴耦合），并且统一了：

- 标准 TEMPO（对角、对易耦合，`ADT` + 部分影响泛函）；
- 对角但非对易的多浴耦合；
- 非对角（共轭对）耦合（`PT` + 平移不变影响泛函）；
- 实时间、虚时间、以及混合（Kadanoff–Baym）轮廓上的演化；
- 含时系统-浴耦合。

---

> 全部导出符号的逐项 API 参考见 [api_reference.md](api_reference.md)。

## 目录

1. [安装](#1-安装)
2. [核心概念](#2-核心概念)
3. [快速上手](#3-快速上手)
4. [核心组件与 API](#4-核心组件与-api)
5. [超参数与误差来源](#5-超参数与误差来源)
6. [代码结构](#6-代码结构)
7. [与文献的对应关系](#7-与文献的对应关系)

---

## 1. 安装

本工具包依赖以下包（见 `Project.toml`）：

| 包 | 作用 |
|---|---|
| `ImpurityModelBase` | 定义谱密度、浴（`bosonicbath`）、玻色算符等基础对象 |
| `QuAPI` | 提供 `ContourIndex`、`branch`、`index` 等轮廓基础类型 |
| `TensorOperations` | 张量网络缩并 |
| `Polynomials`、`LinearAlgebra`、`Statistics`、`TupleTools`、`Logging` | 通用数值与线性代数工具 |

在 Julia 中激活项目后即可使用：

```julia
using TEMPO, ImpurityModelBase, LinearAlgebra
```

> 注意：`spectrum`、`bosonicbath`、`bosonaoperator`、`bosondensityoperator`、`Leggett` 等函数定义在 `ImpurityModelBase` 中，需要一并 `using`。

---

## 2. 核心概念

### 2.1 量子杂质问题（QIP）

考虑一个"杂质"系统 $\hat{H}_S$ 线性耦合到无相互作用的玻色浴：

```math
\hat{H} = \hat{H}_S + \hat{H}_{\text{int}},
\qquad
\hat{H}_{\text{int}} = \hat{H}_{\text{hyb}} + \hat{H}_B .
```

- **对角（diagonal）耦合**（原始 TEMPO/QuAPI）：

  ```math
  \hat{H}_{\text{hyb}} = \sum_{l,k} \hat{A}_l\, (V_{l,k} \hat{b}^\dagger_{l,k} + \mathrm{H.c.}),
  ```
  其中 $\hat{A}_l$ 是厄米算符，耦合项中 $\hat{A}_l$ 只与 $\hat{b}^\dagger+\hat{b}$ 的线性组合成对出现。

- **非对角（off-diagonal）耦合**（本文献的推广）：

  ```math
  \hat{H}_{\text{hyb}} = \sum_{l,k} \left( V_{l,k} \hat{A}_l \hat{b}^\dagger_{l,k} + \mathrm{H.c.} \right),
  ```
  其中 $\hat{A}_l$ 可以是非厄米算符（例如 Jaynes–Cummings 型 $\hat{A}=\hat{\sigma}_-$）。

非对角耦合无法通过重新组合化为"对角且非对易"的情形，需要新的框架。

### 2.2 Feynman–Vernon 影响泛函（IF）

TEMPO 类方法的关键出发点是对浴求迹后得到的 Feynman–Vernon 影响泛函。对非对角耦合，它在 Keldysh 轮廓 $C$ 上具有算符路径形式：

```math
\mathcal{I}[\hat{A}^\dagger, \hat{A}]
= \mathcal{T}_C \exp\left[ -\int_C \mathrm{d}\tau \int_C \mathrm{d}\tau'\,
    \hat{A}^\dagger(\tau)\, \Delta(\tau,\tau')\, \hat{A}(\tau') \right],
```

其中**混合化函数（hybridization function）**由谱密度给出：

```math
\Delta(\tau,\tau') = i \int \mathrm{d}\omega\, J(\omega)\, D_\omega(\tau,\tau'),
\qquad
J(\omega) = \sum_k |V_k|^2 \delta(\omega - \omega_k).
```

对角耦合时该 IF 对应经典配分函数（可表示为 MPS / ADT）；非对角耦合时它对应一个有效量子多体哈密顿量的热态 $\mathrm{e}^{-\hat{H}_{\text{eff}}}$，需要表示为 **MPO（即 PT）**。

### 2.3 两种张量网络：ADT 与 PT

| 对象 | 全称 | 表示 | 适用耦合 |
|---|---|---|---|
| `ADT` | 增强密度张量（Augmented Density Tensor） | MPS | 对角耦合（原始 TEMPO） |
| `ProcessTensor` (`PT`) | 过程张量 | MPO | 对角非对易、非对角耦合（文献扩展） |

一个 PT 可通过对相邻位点施加 3D copy 张量系统地转换为 ADT（文献 Fig. 4）。

### 2.4 轮廓（Contour）

本工具包支持三种时间轮廓，通过 `contour` 关键字选择：

- `contour=:real`（等价于 `:Keldysh`）：实时间演化，系统初态为 $\hat{\rho}_S \otimes \hat{\rho}_B$；
- `contour=:imag`：虚时间演化（$0\to\beta$），对应有限温度配分函数与 Matsubara 关联函数；
- `contour=:mixed`（等价于 `:Kadanoff`）：L 形 Kadanoff–Baym 轮廓（虚时间 + 实时间混合）。

---

## 3. 快速上手

### 3.1 标准自旋玻色子模型（对角耦合，Rabi 型，ADT）

对应 `tutorials/spinboson/rabitype.jl` 与 `benchmark/sb.jl`。系统是二能级自旋，$\hat{A}=\hat{\sigma}_x/2$（对角耦合），浴谱取 sub-Ohmic 谱。

```julia
using TEMPO, ImpurityModelBase, LinearAlgebra

# 截断方案：最大键维 chi，奇异值阈值 ϵ
trunc = truncdimcutoff(D=chi, ϵ=1.0e-12, add_back=0)

# 1. 定义实时间轮廓格点
lattice = ADTLattice(N=Nt, δt=δt, contour=:real)

# 2. 定义系统算符（σ_x 用于系统哈密顿量；σ_z 的对角元用于对角耦合算符 A）
x = Matrix{ComplexF64}([0 1; 1 0])
z = Matrix{ComplexF64}([-1 0; 0 1])
zdiag = [z[i,i] for i in 1:size(z,1)]

# 3. 定义浴与混合化样式
hyb  = AdditiveHyb(zdiag)                       # 对角耦合
spec = spectrum(w -> 2π*α * w^s * wc^(1-s), lb=0, ub=wc)   # sub-Ohmic 谱
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

# 4. 构建影响泛函（IF），得到一个 ADT（MPS）
mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)

# 5. 构建裸系统动力学并施加边界条件（初态）
model = ImpurityHamiltonian(Δ .* x)             # 系统哈密顿量 H_S
mpsK  = sysdynamics(lattice, model, trunc=trunc)
mpsK  = boundarycondition!(mpsK, lattice, ρ₀=ρimp)

# 6. 合并得到完整的 ADT（等价于路径积分积分核）
# mps = mult!(mpsK, mpsI, trunc=trunc)
cache = environments(lattice, mpsK, mpsI)       # 预计算环境，用于观测量

# 7. 计算局部观测量（如 ⟨σ_z(t)⟩）
obs = ComplexF64[]
for i in 1:Nt
    pos = index(lattice, i, branch=:+)
    m   = ADTTerm(pos, zdiag)
    v   = expectationvalue(m, cache)            # 已除以 Z，得到归一化期望值
    push!(obs, v)
end
```

### 3.2 JC 型非对角耦合（PT 框架）

对应 `tutorials/spinboson/jctype.jl`。系统通过 $\hat{A}=\hat{\sigma}_-/2$（共轭对 $\hat{A}^\dagger,\hat{A}$）耦合到浴，这是文献中的**非对角耦合**情形，必须使用 PT 框架与平移不变影响泛函算法。

```julia
lattice = PTLattice(N=Nt, δt=δt, contour=:real)   # 注意是 PTLattice

hyb  = NonDiagonalHyb(sp)                        # 非对角耦合：op*a + op'*a'
spec = spectrum(w -> 2π*α * w^s * wc^(1-s), lb=0, ub=wc)
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

# IF 构建算法：平移不变 IF（XTRG 式），
# 使用 DMRG 型 MPO-MPO 乘法 + Prony 指数展开
algmult  = DMRGMult1(trunc, maxiter=10)
algexpan = PronyExpansion(n=n, tol=1.0e-8, verbosity=2)
alg = TranslationInvariantIF(k=k, fast=true, algmult=algmult, algexpan=algexpan, verbosity=2)

mpsI = hybriddynamics(lattice, corr, hyb, alg)   # 得到 ProcessTensor (MPO)

# 系统动力学 + 初态
model = ImpurityHamiltonian(Δ .* z)
mpsK  = sysdynamics(lattice, model, trunc=trunc)
mps   = mult!(mpsK, mpsI, trunc=trunc)

# 观测量：PT 实时间需要把初态密度矩阵传入环境
cache = environments(lattice, mps, ρ₀=ρimp)

obs = ComplexF64[]
for i in 1:Nt
    ind = ContourIndex(i, :+)
    m   = ContourOperator(ind, x)               # 任意系统算符 x
    v   = expectationvalue(m, cache)
    push!(obs, v)
end
```

### 3.3 玻色杂质：虚时间演化与 Matsubara 关联函数

对应 `benchmark/independentbosons.jl` 与 `bosonicimpurity.jl`。此时杂质本身是玻色模（局域 Hilbert 空间截断为 `d`），在虚时间轮廓上计算 Matsubara Green 函数 $\langle \mathcal{T}_\tau \hat{a}(\tau)\hat{a}^\dagger(0)\rangle$。

```julia
lattice = ADTLattice(N=N, δτ=δτ, contour=:imag)      # 虚时间轮廓
# 或对非对角耦合：PTLattice(...)

a = bosonaoperator(d=d)                              # 玻色湮灭算符（ImpurityModelBase）
n = bosondensityoperator(d=d)

hyb  = AdditiveHyb(diag(n))                          # 对角耦合
spec = Leggett(d=1, ωc=1)                            # 或自定义 spectrum(...)
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

mpsI  = hybriddynamics(lattice, corr, hyb, trunc=trunc)
model = ImpurityHamiltonian(ϵ_d .* n)                # H_S
mpsK  = sysdynamics(lattice, model, trunc=trunc)
Zval  = integrate(mpsK, mpsI)                        # 配分函数 Z

# 两点关联函数：把算符插入到系统动力学中
op1, op2 = [0 0; 1 0], [0 1; 0 0]
c1, c2   = ContourIndex(1), ContourIndex(1)
ct       = ContourOperator([c1, c2], [op1, op2])
mpsK2    = sysdynamics(lattice, model, ct, trunc=trunc)
v        = integrate(mpsK2, mpsI) / Zval
```

对**非对角耦合**的玻色杂质（如 `bosonicimpurity.jl`）：

```julia
lattice = PTLattice(N=N, δt=δt, d=d, contour=:real)
hyb  = NonDiagonalHyb(a')
alg  = TranslationInvariantIF(k=5, fast=true,
                             algmult=DMRGMult1(trunc, initguess=:rand),
                             algexpan=PronyExpansion(n=20, tol=1.0e-8))
mpsI = hybriddynamics(lattice, corr, hyb, alg)
```

### 3.4 含时系统-浴耦合

对含时耦合（见 `src/tdinfluencefunctional/`），使用含时版本的混合化样式：`AdditiveTdHyb`、`NonAdditiveTdHyb`、`NonDiagonalTdHyb`。它们接受一个函数 `f(t)` 描述耦合强度的时间依赖：

```julia
hyb = NonDiagonalTdHyb(op, t -> f(t))   # 或 AdditiveTdHyb / NonAdditiveTdHyb
```

---

## 4. 核心组件与 API

### 4.1 截断方案 `TruncationScheme`

张量网络计算中压缩键维的截断方案（`src/auxiliary/truncation.jl`）：

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

### 4.2 格点 `ADTLattice` / `PTLattice`

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

### 4.3 系统算符 `AbstractImpurityOperator`

| 类型 | 构造 | 说明 |
|---|---|---|
| `ImpurityHamiltonian` | `ImpurityHamiltonian(m::Matrix)` | 幺正演化 $e^{\mp i \hat{H}_S \delta t}$ |
| `ImpurityLindbladian` | `ImpurityLindbladian(H, jumpops)` | 含耗散（Lindblad）演化，`jumpops` 为跃迁算符列表 |

由 `sysdynamics(lattice, model, trunc=...)` 生成裸系统动力学张量 `mpsK`（ADT 或 PT）。

### 4.4 混合化样式 `HybridizationStyle`

描述系统-浴耦合的具体形式（`src/influencefunctional/influencefunctional.jl`、`src/tdinfluencefunctional/`）：

| 类型 | 耦合形式 | 约束 | 框架 |
|---|---|---|---|
| `AdditiveHyb(op)` | $\hat{A}(\hat{b}^\dagger+\hat{b})$ 型，对角 | `op` 为向量或对角矩阵（厄米） | ADT |
| `NonAdditiveHyb(op)` | $\hat{A}(\hat{a}+\hat{a}^\dagger)$，$\hat{A}=\hat{A}^\dagger$ | `op` 厄米矩阵 | PT |
| `NonDiagonalHyb(op)` | $\hat{A}\hat{a} + \hat{A}^\dagger \hat{a}^\dagger$ | `op` 任意方阵 | PT（文献的核心情形） |

含时版本（`TdHybridizationStyle`）：`AdditiveTdHyb`、`NonAdditiveTdHyb`、`NonDiagonalTdHyb`。

> 对于对角耦合，`hybriddynamics` 既可用 `PartialIF`（逐位点乘，每因子键维为 2，见文献 [SciPost Phys. Core 7, 063 (2024)]），也可用 `TranslationInvariantIF`。对非对角耦合只能使用 `TranslationInvariantIF`。

### 4.5 影响泛函构建算法 `InfluenceFunctionalAlgorithm`

```julia
# 1) 部分影响泛函（Partial IF）：仅适用于对角耦合 AdditiveHyb
alg = PartialIF(trunc=trunc)

# 2) 平移不变影响泛函（TTI-IF，XTRG 式），文献推荐
alg = TranslationInvariantIF(;
    algexpan = PronyExpansion(n=20, tol=1.0e-8, verbosity=2),  # 混合化函数指数展开
    algevo   = WII(),        # 或 WI()、ComplexStepper()、FirstOrderStepper()
    algmult  = DMRGMult1(trunc, initguess=:rand, maxiter=10),   # 或 SVDCompression(trunc)
    k        = 7,            # XTRG 步数：时间步长 1/2^k
    fast     = true,         # true：树形二分方案（约 k 次乘法）；false：串行 2^k-1 次
    verbosity= 0,
)
```

- `algexpan`：`ExponentialExpansionAlgorithm`，包括 `PronyExpansion`、`DeterminedPronyExpansion`（`exponential_expansion`、`expansion_error` 可用于误差分析）；
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

### 4.6 浴与关联函数（`ImpurityModelBase`）

```julia
spec = spectrum(f, lb=0, ub=wc)          # 由函数定义谱密度 J(ω)
spec = Leggett(d=1, ωc=5, α=0.1)         # 预定义 Leggett 谱
bath = bosonicbath(spec, β=β)            # 玻色浴（有限温度 β）
corr = correlationfunction(bath, lattice) # 离散化到格点上的关联函数 Δ(τ,τ')
```

`correlationfunction` 会根据格点轮廓自动调用 `Δt`（实时间）、`Δτ`（虚时间）或 `Δm`（混合轮廓）。

### 4.7 ADT / PT 张量操作

- `mult(a, b, algmult)` / `mult!(a, b, trunc)`：张量网络乘法（元素级乘积 / MPO 乘法）；
- `integrate(lattice, args...)` / `integrate(mpsA, mpsB)`：求配分函数（路径积分求和）；
- `apply!(term, mps)`：应用局域算符；
- `canonicalize!`、`leftorth!`、`rightorth!`：正交化（可指定 `Orthogonalize`）；
- `bond_dimension(mps)`、`bond_dimensions(mps)`：键维查询；
- `distance(mps1, mps2)` / `distance2`：两个张量的距离（相对误差验证用）；
- `randomadt` / `randompt`：随机张量（测试用）；`vacuumstate(lattice)`：真空态。

### 4.8 观测量

**实时间（初始密度矩阵 $\hat\rho_0$）**

```julia
# PT 框架（非对角耦合）
cache = environments(lattice, mps, ρ₀=ρimp)
v = expectationvalue(ContourOperator(ContourIndex(i, :+), op), cache)

# ADT 框架（对角耦合）
cache = environments(lattice, mpsK, mpsI)
v = expectationvalue(ADTTerm(index(lattice, i, branch=:+), zdiag), cache)
```

**虚时间 / 混合轮廓**

```julia
cache = environments(lattice, mpsK, mpsI)   # ADT
cache = environments(lattice, mps)          # PT
v = expectationvalue(ContourOperator(...), cache)
```

辅助函数：`Zvalue(cache)`（配分函数）、`Zvalue2(cache)`、`TransferMatrix`（转移矩阵）、
`correlation(lattice, model, op, mpsI[, ρ0])`（两点关联函数）、`heatcurrents`。

### 4.9 MPO 哈密顿量（长程相互作用工具）

`src/mpohamiltonian/` 提供了构造长程（指数衰减）相互作用哈密顿量 MPO 的独立工具，用于其他一维量子多体问题：

```julia
# SchurMPOTensor：把 [局域项 + 指数衰减长程项] 编码为紧凑 MPO 站点张量
h = SchurMPOTensor(h1, h2s)    # h2s 为 ExponentialDecayTerm / GenericDecayTerm / PowerlawDecayTerm 列表
mpo = MPOHamiltonian([h for _ in 1:L])
tensors = tompotensors(mpo)              # 转成稠密 MPO 站点张量
tensors2 = timeevompo(tensors, dt, WII())   # 时间演化（WI / WII / ComplexStepper / FirstOrderStepper）
```

---

## 5. 超参数与误差来源

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

---

## 6. 代码结构

```
src/
├── TEMPO.jl / includes.jl      # 模块定义与导出符号
├── auxiliary/                  # 截断、DMRG 乘法、张量操作、正交化等基础工具
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

tutorials/spinboson/            # 教程脚本：rabitype（标准 TEMPO）与 jctype（非对角耦合）
benchmark/                      # 基准算例：自旋玻色子、独立玻色模、实时间等
test/                           # 测试套件（含与 Exact Diagonalization 的对照）
```

---

## 7. 与文献的对应关系

| 文献内容 | 工具包实现 |
|---|---|
| 对角 + 对易耦合（标准 TEMPO） | `ADTLattice` + `AdditiveHyb` + `PartialIF`/`TranslationInvariantIF` |
| 对角 + 非对易耦合（多浴） | 多个 `AdditiveHyb` 的 IF 相乘（`mult!`） |
| 非对角耦合（共轭对，核心推广） | `PTLattice` + `NonDiagonalHyb` + `TranslationInvariantIF` |
| 实时间 Keldysh 轮廓 | `contour=:real` |
| 虚时间轮廓 | `contour=:imag` |
| L 形 Kadanoff–Baym 轮廓 | `contour=:mixed`（`MixedPTLattice`/`MixedADTLattice`） |
| QuAPI 离散化（附录 C） | `correlationfunction(bath, lattice)` 中的 `Δt`/`Δτ`/`Δm` |
| 混合化函数的指数展开（附录 D） | `PronyExpansion`/`DeterminedPronyExpansion`（`algexpan`） |
| XTRG 构造有效热态（文献 Fig. 3） | `TranslationInvariantIF(k=..., fast=...)` |
| PT 中系统哈密顿量的吸收（文献 Fig. 2b） | `sysdynamics(lattice, model)` + `mult!` |
| 观测量计算（文献 Fig. 1d, 1e） | `environments` + `expectationvalue` |
| JC 自旋玻色子模型 | `tutorials/spinboson/jctype.jl` |
| 标准自旋玻色子模型 | `tutorials/spinboson/rabitype.jl`、`benchmark/sb.jl` |
| 非相互作用/相互作用玻色杂质 | `benchmark/independentbosons.jl`、`bosonicimpurity.jl` |

---

*本文档基于 `docs/GuoChen2026a.pdf`（Phys. Rev. B 114, 125413 (2026)）整理。*
