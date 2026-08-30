# TEMPO.jl API 参考（导出符号）

本文档逐一说明 `src/TEMPO.jl` 中 `export` 的全部符号（类型、函数、常量），按功能分类组织。所有签名均与源码核对。

使用方式：

```julia
using TEMPO, ImpurityModelBase, LinearAlgebra
```

> 说明：
> - `scalartype`、`space_l`、`space_r` 由 `TensorOperations` 提供，经 TEMPO 重导出；
> - `branch`、`index` 由 `QuAPI` 提供，经 TEMPO 重导出；
> - `phydim`、`spectrum`、`bosonicbath` 等在 `ImpurityModelBase` 中定义（不在此清单内）。

---

## 目录

- [1. 截断方案 TruncationScheme](#1-截断方案)
- [2. 指数展开 Exponential Expansion](#2-指数展开)
- [3. 张量分解与矩阵操作](#3-张量分解与矩阵操作)
- [4. 轮廓索引 ContourIndex](#4-轮廓索引)
- [5. MPO 哈密顿量](#5-mpo-哈密顿量)
- [6. 增强密度张量 ADT](#6-增强密度张量-adt)
- [7. 过程张量 ProcessTensor](#7-过程张量-processtensor)
- [8. 局域算符项](#8-局域算符项)
- [9. 格点 Lattice](#9-格点-lattice)
- [10. PT 格点求值函数](#10-pt-格点求值函数)
- [11. 影响泛函 IF](#11-影响泛函-if)
- [12. 边界条件与初态](#12-边界条件与初态)
- [13. 系统模型 Models](#13-系统模型-models)
- [14. 观测量 Observables](#14-观测量-observables)

---

## 1. 截断方案

```julia
abstract type TruncationScheme end
```

所有截断方案的抽象基类型，被 `trunc` 关键字接受。

### NoTruncation

```julia
struct NoTruncation <: TruncationScheme end
NoTruncation()
```

不截断（保留全部奇异值）。

### truncdimcutoff — 推荐的截断方案

```julia
truncdimcutoff(D=χ, ϵ=ε, add_back=0)  # → TruncationDimCutoff
struct TruncationDimCutoff <: TruncationScheme
    D::Int          # 最大键维
    ϵ::Float64      # 奇异值截断阈值
    add_back::Int   # 是否额外保留若干维
end
```

同时限制最大键维 `D` 与奇异值阈值 `ϵ`，推荐用于所有计算。

### truncdim

```julia
truncdim(D)     # → TruncateDim
struct TruncateDim <: TruncationScheme
    D::Int
end
```

只按最大键维截断。

### trunccutoff

```julia
trunccutoff(ϵ)  # → TruncateCutoff
struct TruncateCutoff <: TruncationScheme
    ϵ::Float64
end
```

只按奇异值阈值截断。

### renyi_entropy

```julia
renyi_entropy(v::AbstractVector{<:Real}; α::Real=1)
```

计算归一化向量 `v`（满足 `sum(v)==1`）的 Rényi 熵。`α=1` 退化为 von Neumann 熵。

---

## 2. 指数展开

用于将混合化函数等数据 $f(x)$（$1\le x\le N$）展开为指数和 $\sum_i \alpha_i \beta_i^{x}$。

### OverDeterminedProny / DeterminedProny

```julia
struct OverDeterminedProny <: AbstractPronyExpansion
    n::Int; tol::Float64; verbosity::Int; stepsize::Union{Int, Nothing}
end
OverDeterminedProny(; n::Int=10, tol::Real=1.0e-8, verbosity::Int=1, stepsize::Union{Int,Nothing}=1)

struct DeterminedProny <: AbstractPronyExpansion
    n::Int; tol::Float64; verbosity::Int; stepsize::Union{Int, Nothing}
end
DeterminedProny(; n::Int=10, tol::Real=1.0e-8, verbosity::Int=1, stepsize::Union{Int,Nothing}=1)
```

- `n`：最大展开项数（迭代上限）；
- `tol`：相对误差阈值（`tol*norm(f)`），达到即提前收敛；
- `verbosity`：输出级别（≥1 打印不收敛警告，≥2 打印收敛信息）；
- `stepsize`：均匀采样步长（默认 1，即在 `f[1:s:end]` 上拟合并换算回原采样）；`nothing` 时自动选步长。

两者区别：`DeterminedProny` 使用确定型 Prony 方法（无最小二乘），`OverDeterminedProny` 使用最小二乘 Prony 方法。

### exponential_expansion

```julia
exponential_expansion(f::Vector{<:Number}, alg::ExponentialExpansionAlgorithm)
exponential_expansion(f::Vector{<:Number}; alg::ExponentialExpansionAlgorithm=OverDeterminedProny())
exponential_expansion(f, L::Int, alg::ExponentialExpansionAlgorithm)
exponential_expansion(f, L::Int; alg::ExponentialExpansionAlgorithm=OverDeterminedProny())
```

返回展开系数与底数 `(αs::Vector, βs::Vector)`，满足 `f(x) ≈ Σᵢ αᵢ * βᵢˣ`（`x = 1,2,…,N`）。后两种形式接受函数 `f(k)` 并自动采样 `k=1:L`。

### expansion_error

```julia
expansion_error(f::Vector{<:Number}, p::Vector{<:Number})
expansion_error(f::Vector{<:Number}, coeffs::Vector{<:Number}, alphas::Vector{<:Number})
```

计算指数展开在采样点上的 2-范数误差，用于判断展开质量。

---

## 3. 张量分解与矩阵操作

### OrthogonalFactorizationAlgorithm

```julia
abstract type OrthogonalFactorizationAlgorithm <: FactorizationAlgorithm end
# 具体算法：QR(), QRpos(), LQ(), LQpos(), SVD(), SDD(), Polar()
```

正交分解算法类型，用于 `leftorth!` / `rightorth!`。

### leftorth! / rightorth!

```julia
leftorth!(A::StridedMatrix; alg=QRpos(), atol=0.0)
leftorth!(A::AbstractArray{T,N}, left::NTuple, right::NTuple; alg=QRpos(), atol=0.0)

rightorth!(A::StridedMatrix; alg=LQpos(), atol=0.0)
rightorth!(A::AbstractArray{T,N}, left::NTuple, right::NTuple; alg=LQpos(), atol=0.0)
```

- 矩阵版本返回 `(Q, R)` / `(L, Q)`；
- 张量版本将 `left` 维合并为左因子，`right` 维合并为右因子，返回带重排形状的分解结果。

### permute

```julia
permute(m::AbstractArray, perm)            # → PermutedDimsArray
permute(m::AbstractArray, left, right)
```

不拷贝数据的维序重排（视图）。`permute(m, left, right)` 等价于 `permute(m, (left..., right...))`。

### tsvd!

```julia
tsvd!(a::StridedMatrix, workspace=similar(a, length(a)); trunc::TruncationScheme=NoTruncation())
tsvd!(a::StridedArray{T,N}, left::NTuple, right::NTuple, workspace=similar(a, length(a)); trunc::TruncationScheme=NoTruncation())
```

带截断的紧凑 SVD。返回 `(U, S, V, err)`，其中 `err` 为截断误差。内部先用 divide-and-conquer（`gesdd`），失败时回退到 `gesvd`（数值稳定性）。

### scalartype

```julia
scalartype(x)
```

返回对象（张量/类型）的元素标量类型。

### space_l / space_r / bond_dimension / bond_dimensions / phydim / phydims / scaling

```julia
space_l(x, i)      # 第 i 维"左"辅助空间维数（第 1 指标）
space_r(x, i)      # 第 i 维"右"辅助空间维数（第 3 指标）
bond_dimension(mps)     # 最大键维 max(bond_dimensions)
bond_dimensions(mps)    # 各键的键维向量
phydim(x[, i])          # 物理维数 d（i 省略时对均匀张量/格点）
phydims(mps)            # 各站点物理维数
scaling(mps)            # 张量整体缩放因子（用于数值稳定）
```

---

## 4. 轮廓索引

```julia
struct ContourIndex
    j::Int
    branch::Symbol    # ∈ (:+, :-, :τ)
end
ContourIndex(j::Int, branch::Symbol)
ContourIndex(j::Int; branch::Symbol=:τ)
```

轮廓索引，`branch` 含义：`:+` 前向实时间分支，`:-` 后向实时间分支，`:τ` 虚时间分支。

```julia
branch(x::ContourIndex)   # 返回分支符号
index(lattice, i, branch) # (时间步, 分支) → 格点位点位置
branches(lattice)         # 格点包含的分支元组，如 (:+, :-) / (:τ,) / (:+, :-, :τ)
```

`ContourIndex` 定义了按轮廓顺序的偏序 `<`（Keldysh 次序）。

---

## 5. MPO 哈密顿量

一维长程相互作用哈密顿量的稀疏 MPO 工具（也用于 IF 中有效哈密顿量的指数化）。

### MPOHamiltonian

```julia
struct MPOHamiltonian{M <: AbstractSparseMPOTensor}
    data::Vector{M}
end
MPOHamiltonian(data::Vector{M})                  # M <: AbstractSparseMPOTensor
MPOHamiltonian(data::Vector{Matrix{Any}})        # 由 SparseMPOTensor 矩阵构造
```

站点张量为稀疏形式的哈密顿量，支持 `length` / `getindex(h, i)` / `getindex(h, i, j, k)` 等。

### tompotensors

```julia
tompotensors(h::MPOHamiltonian)
tompotensors(h::MPOHamiltonian{<:SparseMPOTensor}; rowl::Int=1, colr::Int=1)
```

将 `MPOHamiltonian` 展开为稠密 MPO 站点张量数组（`Vector{Array{T,4}}`，即 `ProcessTensor` 可接受的数据）。

### WI / WII / FirstOrderStepper / ComplexStepper

```julia
abstract type TimeEvoMPOAlgorithm <: MPSAlgorithm end
abstract type FirstOrderStepper <: TimeEvoMPOAlgorithm end

WI(; tol::Real=Defaults.tol, maxiter::Int=Defaults.maxiter)
WII(; tol::Real=Defaults.tol, maxiter::Int=Defaults.maxiter)

struct ComplexStepper{F<:FirstOrderStepper} <: SecondOrderStepper
    stepper::F
end
```

时间演化步进器（"WI / WII" 命名取自 [Zaletel et al., arXiv:1407.1832]）：
- `WI`：一阶，局域项用 $1+\mathrm{d}t\,D$ 近似；
- `WII`：一阶，局域项精确指数化 $e^{\mathrm{d}t\,D}$（精度更高）；
- `ComplexStepper(stepper)`：二阶（Suzuki 分解），对复数时间步长 `dt` 做了特殊处理。

### timeevompo

```julia
timeevompo(mpo, dt::Number, alg::TimeEvoMPOAlgorithm)
```

把 (平移不变) MPO 演化一步 `dt`：
- 对 `SchurMPOTensor` 输入返回新的 `SchurMPOTensor`；
- 对 `FirstOrderStepper`（`WI`/`WII`）返回单个 MPO；
- 对 `ComplexStepper` 返回一对 MPO（`mpo⁺, mpo⁻`，对应 $e^{-iH\,\mathrm{d}t/2}$ 与 $e^{+iH\,\mathrm{d}t/2}$）。

### SchurMPOTensor

```julia
struct SchurMPOTensor <: AbstractSparseMPOTensor end

SchurMPOTensor(h1::AbstractMatrix, h2s::Vector)   # h1 局域项, h2s 长程项
SchurMPOTensor(h2s::Vector)                        # 无局域项，等价 h1=0
```

把 [局域项 + 指数衰减长程项] 编码为 Schur 形式的稀疏站点张量（外层矩阵尺寸 $(N+2)\times(N+2)$，`N=length(h2s)`）。

### ExponentialDecayTerm

```julia
struct ExponentialDecayTerm <: AbstractLongRangeTerm
    a::M1; m::M2; b::M3; α::T; coeff::T
end
ExponentialDecayTerm(a, b; middle=_eye(size(a,1)), α::Number=1., coeff::Number=1.)
```

指数衰减长程项 $\text{coeff}\cdot \hat{a}\otimes(\alpha\hat{m})^{\otimes n}\otimes \hat{b}$（`m` 为中间传播算符，`α` 为衰减比）。

### GenericDecayTerm / PowerlawDecayTerm

```julia
GenericDecayTerm(a, b; f, middle=_eye(size(a,1)), coeff::Number=1.)
# f: 衰减函数，作用于中间维

PowerlawDecayTerm(a, b; α::Number=1., kwargs...)   # 等价 GenericDecayTerm(f=x->x^α)
```

一般衰减函数的长程项（`PowerlawDecayTerm` 为幂律衰减 $n^{-\alpha}$）。

---

## 6. 增强密度张量 ADT

### ADT

```julia
struct ADT{T<:Number, R<:Real} <: Dense1DTN{T}
    data::Vector{Array{T,3}}
    s::Vector{Union{Missing, Vector{R}}}   # 奇异向量（Schmidt 谱）
    scaling::Ref{Float64}
end

ADT(data::Vector{<:DenseMPSTensor}; scaling::Real=1)
ADT(::Type{T}, ds::AbstractVector{Int})        # 均匀全 1 张量，物理维 ds
ADT(::Type{T}, L::Int; d::Int=2)
ADT(ds::AbstractVector{Int}); ADT(L::Int; d::Int=2)
```

增强密度张量，站点张量秩 3（MPS）。边界键维固定为 1。常用属性：`length`、`getindex/setindex!`、`bond_dimension`、`phydim`、`scaling`、`copy`、`complex`。

### randomadt

```julia
randomadt(::Type{T}, ds::AbstractVector{Int}; D::Int)
randomadt(::Type{T}, L::Int; D::Int, d::Int=2)
randomadt(ds; kwargs...); randomadt(L; kwargs...)
```

生成随机 ADT（键维 `D`）。

### 正交性检查

```julia
isleftcanonical(a::ADT; kwargs...)
isrightcanonical(a::ADT; kwargs...)
iscanonical(a::ADT; kwargs...)      # 右正交且奇异向量为真实 Schmidt 数
```

### distance / distance2

```julia
distance(x, y)   # = sqrt(distance2(x, y))，两个张量/向量的 2-范数距离
distance2(x, y)  # = ‖x‖² + ‖y‖² - 2 Re⟨x,y⟩（取绝对值）
```

用于验证乘法精度（如 `hybriddynamics` 与朴素实现对比）。

### Orthogonalize

```julia
struct Orthogonalize{A<:Union{QR,SVD}, T<:TruncationScheme} <: MatrixProductOrthogonalAlgorithm
    orth::A; trunc::T; normalize::Bool; verbosity::Int
end
Orthogonalize(; alg::Union{QR,SVD}=SVD(), trunc::TruncationScheme=NoTruncation(),
              normalize::Bool=false, verbosity::Int=0)
```

正交化算法配置：`alg` 选 QR 或 SVD（SVD 支持截断），`normalize` 是否归一化整体范数。

### leftorth! / rightorth! / canonicalize!

```julia
leftorth!(psi::ADT; alg::Orthogonalize=Orthogonalize())
rightorth!(psi::ADT; alg::Orthogonalize=Orthogonalize())
canonicalize!(psi::ADT; alg::Orthogonalize=Orthogonalize(trunc=DefaultTruncation, normalize=false))
canonicalize(psi::ADT; kwargs...)   # 非就地版本（deepcopy）
```

- `leftorth!`：变为左正则形式；
- `rightorth!`：变为右正则形式；
- `canonicalize!`：先无截断左正交，再按 `alg` 右正交并截断，返回规范形式并更新奇异向量 `psi.s`。

### mult / mult!

```julia
mult(x::ADT, y::ADT, alg::DMRGAlgorithm)       # 元素级（Hadamard 型）乘积，返回新对象
mult(x::ADT, y::ADT, alg::SVDCompression)
mult!(x::ADT, y::ADT, alg::DMRGAlgorithm)      # 就地：x ← x ⊙ y
mult!(x::ADT, y::ADT, alg::SVDCompression)
mult!(x::ADT, y::ADT; trunc::TruncationScheme=DefaultITruncation, verbosity::Int=0)
```

ADT 的乘法（张量元素级乘积），自动调整缩放因子 `scaling(x)*scaling(y)`。`alg` 决定压缩方式：
- `SVDCompression`：直接 SVD 压缩；
- `DMRGMult1`（推荐，默认）：迭代式 DMRG 乘法。

### DMRGMult1

```julia
struct DMRGMult1 <: DMRGMultAlgorithm
    trunc::TruncationDimCutoff
    maxiter::Int; tol::Float64; initguess::Symbol; verbosity::Int; callback::Function
end
DMRGMult1(trunc::TruncationDimCutoff;
          maxiter::Int=5, tol::Float64=1.0e-12, initguess::Symbol=:svd,
          verbosity::Int=0, callback::Function=Returns(nothing))
DMRGMult1(; trunc::TruncationDimCutoff=DefaultITruncation, kwargs...)
```

单站点迭代 DMRG 乘法算法：
- `maxiter`：最大 sweep 次数；`tol`：相对误差收敛阈值；
- `initguess ∈ (:svd, :pre, :rand)`：初始猜测（SVD 初始、预置 x、随机）；
- 返回结果的 `scaling` 为两输入 scaling 之积。

---

## 7. 过程张量 ProcessTensor

### ProcessTensor

```julia
struct ProcessTensor{T<:Number, R<:Real} <: Dense1DTN{T}
    data::Vector{Array{T,4}}
    s::Vector{Union{Missing, Vector{R}}}
    scaling::Ref{Float64}
end
ProcessTensor(data::Vector{<:DenseMPOTensor}; scaling::Real=1)
ProcessTensor(::Type{T}, ds::AbstractVector{Int})
ProcessTensor(::Type{T}, L::Int; d::Int=2)
ProcessTensor(ds); ProcessTensor(L; d=2)
```

过程张量（MPO），站点张量秩 4（左辅助、输入物理、右辅助、输出物理）。用于非对角耦合。支持 `*`（元素级乘积）、`+`/`-`（MPO 相加，右边界需维数匹配）、`copy`、`complex`、`bond_dimension`、`phydim`、`scaling`。

### randompt

```julia
randompt(::Type{T}, ds::Vector{Int}; D::Int)
randompt(::Type{T}, L::Int; d::Int=2, D::Int)
randompt(ds; kwargs...); randompt(L; kwargs...)
```

生成随机 ProcessTensor（辅助维 `D`）。

### 正交化（与 ADT 相同）

```julia
leftorth!(h::ProcessTensor; alg::Orthogonalize=Orthogonalize())
rightorth!(h::ProcessTensor; alg::Orthogonalize=Orthogonalize())
canonicalize!(h::ProcessTensor; alg::Orthogonalize=Orthogonalize(trunc=DefaultTruncation))
isleftcanonical / isrightcanonical / iscanonical(h::ProcessTensor)
```

### mult / mult!

```julia
mult(x::ProcessTensor, y::ProcessTensor, alg::DMRGAlgorithm)
mult!(x::ProcessTensor, y::ProcessTensor, alg::DMRGAlgorithm)
mult!(x::ProcessTensor, y::ProcessTensor; trunc=DefaultITruncation, verbosity=0)
```

MPO-MPO 乘法（非对角 IF 的核心运算），压缩算法同上（`SVDCompression` / `DMRGMult1`）。

---

## 8. 局域算符项

### ADTTerm

```julia
struct ADTTerm{N, T<:Number} <: AbstractADTTerm
    data::Vector{Array{T,3}}
    positions::NTuple{N, Int}
end
ADTTerm(positions::NTuple{N, Int}, data::AbstractArray{T,N})
ADTTerm(p::Int, data::AbstractVector{<:Number})          # 对角算符（单点）
ADTTerm(p::Int, data::DenseMPSTensor)
ADTTerm(positions::NTuple{N, Int}, data::NTuple{N, <:AbstractVector})
```

ADT 上的局域算符项（作用在若干位点上的算符积），位置自动排序。

### apply!

```julia
apply!(x::ADTTerm, mps::ADT)                    # 就地应用，返回 mps
apply!(x::GenericFockTerm, mps::ProcessTensor)
apply!(x::ProdFockTerm, mps::ProcessTensor; aheads=true)
apply!(x::ContourOperator, lat::AbstractPTLattice, mps::ProcessTensor; aheads=true)
```

把局域项乘到张量网络指定位点（对 PT 需要额外给出格点用于索引映射）。

### FockTerm 系列

```julia
abstract type AbstractFockTerm{T<:Number} end

FockTermS(positions::NTuple{N,Int}, data)     # 固定长度 MPO 局域项（秩 4 站点张量）
FockTerm(positions::AbstractVector{Int}, data::AbstractVector{<:AbstractArray{T,4}})
ProdFockTerm(positions::AbstractVector{Int}, data::AbstractVector{<:AbstractMatrix})
ProdFockTerm(pos::Int, data::AbstractMatrix)
```

PT（MPO）上的局域算符项：
- `FockTermS` / `FockTerm`：由站点张量（秩 4）构造的多点 MPO 项；
- `ProdFockTerm`：乘积形式的局域算符项（每个位置一个矩阵），`aheads` 控制插入顺序。

### ContourOperator

```julia
struct ContourOperator{T<:Number}
    indices::Vector{ContourIndex}
    ops::Vector{Matrix{T}}
end
ContourOperator(idx::AbstractVector{ContourIndex}, data::AbstractVector{<:AbstractMatrix})
ContourOperator(p::ContourIndex, data::AbstractMatrix)
```

PT 上的轮廓算符：每个 `ContourIndex` 位置作用一个矩阵（`:-` 分支自动取转置）。用于构造两点关联函数的插入算符。

---

## 9. 格点 Lattice

### 构造入口

```julia
ADTLattice(; contour::Symbol, kwargs...)    # :real/:Keldysh, :imag, :mixed/:Kadanoff
PTLattice(; contour::Symbol, kwargs...)
```

按 `contour` 分发到具体类型：

| 轮廓 | ADT 类型 | PT 类型 |
|---|---|---|
| `:real` / `:Keldysh` | `RealADTLattice`, `RealADTLattice1Order` | `RealPTLattice`, `RealPTLattice1Order` |
| `:imag` | `ImagADTLattice`, `ImagADTLattice1Order` | `ImagPTLattice`, `ImagPTLattice1Order` |
| `:mixed` / `:Kadanoff` | `MixedADTLattice`, `MixedADTLattice1Order` | `MixedPTLattice`, `MixedPTLattice1Order` |

常用关键字（按轮廓而定）：
```julia
ADTLattice(N=Nt, δt=δt, contour=:real)            # 实时间
ADTLattice(N=Nτ, δτ=δτ, contour=:imag)            # 虚时间（d=局域维数）
PTLattice(Nt=Nt, Nτ=Nτ, δt=δt, δτ=δτ, d=d, contour=:mixed)
```

公共接口：`length(lattice)`、`lattice.d`、`lattice.N/Nt/Nτ`、`lattice.δt/δτ`、`lattice.β`、`lattice.t`、`branches(lattice)`、`phydim(lattice)`、`index(lattice, i, branch)`、`lattice[ContourIndex]`、`scalartype`。

### vacuumstate

```julia
vacuumstate(lattice)            # 与格点类型对应的"真空"张量（全 1 均匀态）
vacuumstate(::Type{T}, lattice)
```

ADT 格点返回 `ADT`，PT 格点返回 `ProcessTensor`。

### indexmappings

```julia
indexmappings(lattice)
```

返回轮廓索引 (i, branch) ↔ 位点位置的映射表（具体返回结构随格点类型而异）。

### Fock 排序相关类型

```julia
FockOrdering, ImagFockOrdering, RealFockOrdering, MixedFockOrdering   # 排序方案类型
TimeOrderingStyle; OrderingStyle; LayoutStyle
ImaginaryTimeOrderingStyle; RealTimeOrderingStyle
M2M1; M2m2M1m1; M2M1_m1M1m2M2; TimeLocalLayout                       # 具体排序/布局实现
```

用于控制格点的时间排序与 Fock 空间布局（一般用户无需直接使用，由 `contour` 关键字自动选择）。`OrderingStyle(lattice)`、`LayoutStyle(lattice)`、`ImaginaryTimeOrderingStyle(lattice)`、`RealTimeOrderingStyle(lattice)` 返回对应风格对象。

---

## 10. PT 格点求值函数

### integrate

```julia
integrate(x::ADT)                      # 单张量配分函数
integrate(x::ADT, y::ADT)              # 双张量配分函数
integrate(lat::AbstractPTLattice, x::ProcessTensor[, y::ProcessTensor])
integrate(lattice, mpsK, mpsI)         # 系统 × 影响泛函（常见用法）
```

对张量网络全缩并求和（路径积分求和），返回配分函数 $Z$（标量）。

### rdm

```julia
rdm(lat::RealPTLattice, x::ProcessTensor[, y::ProcessTensor])
```

实时间 PT 的最终约化密度矩阵（输出量子态）。

### quantummap

```julia
quantummap(lat::RealPTLattice, x::ProcessTensor)
```

实时间 PT 的量子映射（多体 Choi 态，未对末端求迹）。

### meanforcestate / mfs

```julia
meanforcestate(lat::ImagPTLattice, x::ProcessTensor[, y::ProcessTensor])
mfs(lat::ImagPTLattice, x::Vararg{ProcessTensor})      # 别名
```

虚时间 PT 的平均力吉布斯态（Mean-Force Gibbs State）。

### correlationfunction（浴关联函数）

```julia
correlationfunction(bath::AbstractBath, lattice::Union{ImagADTLattice1Order, ImagPTLattice1Order})
correlationfunction(bath::AbstractBath, lattice::Union{RealADTLattice1Order, RealPTLattice1Order})
correlationfunction(bath::AbstractBath, lattice::Union{MixedADTLattice1Order, MixedPTLattice1Order})
```

把浴的混合化函数离散化到格点上，返回 `ImagCorrelationFunction` / `RealCorrelationFunction` / `MixedCorrelationFunction`（内部封装 `Δτ` / `Δt` / `Δm`）。`scalartype(corr)` 对实时间返回 `ComplexF64`。

---

## 11. 影响泛函 IF

### 混合化样式

```julia
abstract type HybridizationStyle end
abstract type GeneralHybStyle <: HybridizationStyle end
```

#### AdditiveHyb（对角耦合，标准 TEMPO）

```julia
struct AdditiveHyb <: HybridizationStyle
    op::Vector{Float64}
end
AdditiveHyb(x::AbstractVector{<:Real})
AdditiveHyb(a::AbstractMatrix)     # 要求对角矩阵，取对角元
```

耦合 $\hat{A}(\hat{b}^\dagger+\hat{b})$ 型，`op` 为 $\hat{A}$ 的对角元向量（要求厄米/对易）。

#### NonAdditiveHyb（厄米非对易耦合）

```julia
struct NonAdditiveHyb{T<:Number} <: GeneralHybStyle
    op::Matrix{T}
end
NonAdditiveHyb(a::AbstractMatrix)   # 要求厄米矩阵
NonAdditiveHyb(hyb::AdditiveHyb)
```

耦合 $\hat{A}(\hat{a}+\hat{a}^\dagger)$，$\hat{A}=\hat{A}^\dagger$。

#### NonDiagonalHyb（非对角耦合，文献核心）

```julia
struct NonDiagonalHyb{T<:Number} <: GeneralHybStyle
    op::Matrix{T}
end
NonDiagonalHyb(a::AbstractMatrix)
NonDiagonalHyb(hyb::AdditiveHyb); NonDiagonalHyb(hyb::NonAdditiveHyb)
```

耦合 $\hat{A}\hat{a} + \hat{A}^\dagger\hat{a}^\dagger$，`op` 为任意方阵（可非厄米，如 $\hat{A}=\hat{\sigma}_-$）。

#### pairop

```julia
pairop(b::HybridizationStyle)   # → (op1, op2)
```

返回耦合算符共轭对 $(\hat{A}, \hat{A}^\dagger)$：
- `AdditiveHyb` / `NonAdditiveHyb`：`(op, op)`；
- `NonDiagonalHyb`：`(op, op')`。

### IF 算法

```julia
abstract type InfluenceFunctionalAlgorithm end
```

#### PartialIF

```julia
struct PartialIF <: InfluenceFunctionalAlgorithm
    trunc::TruncationDimCutoff
end
PartialIF(; trunc::TruncationDimCutoff=DefaultITruncation)
```

部分影响泛函：把 IF 写为若干键维 D=2 的 MPS 因子之积（对角耦合专用，见 [SciPost Phys. Core 7, 063 (2024)]）。

#### TranslationInvariantIF

```julia
struct TranslationInvariantIF <: InfluenceFunctionalAlgorithm
    algexpan::ExponentialExpansionAlgorithm
    algevo::TimeEvoMPOAlgorithm
    algmult::DMRGAlgorithm
    k::Int
    fast::Bool
    verbosity::Int
end
TranslationInvariantIF(;
    algexpan::ExponentialExpansionAlgorithm=OverDeterminedProny(n=15, tol=1.0e-4, verbosity=0),
    algevo::TimeEvoMPOAlgorithm=WII(),
    algmult::DMRGAlgorithm=DefaultMultAlg,      # DMRGMult1(DefaultITruncation)
    k::Int=5,
    fast::Bool=true,
    verbosity::Int=0)
```

平移不变 IF（XTRG 式热态构造）：
- `algexpan`：混合化函数指数展开算法；
- `algevo`：有效哈密顿量指数化步进器（`WI`/`WII`/`ComplexStepper`/`FirstOrderStepper`）；
- `algmult`：MPO 乘法压缩算法（`DMRGMult1` 或 `SVDCompression`）；
- `k`：二分步数，最小时间步 $1/2^k$；
- `fast=true` 用树形二分（约 `k` 次乘法），`fast=false` 用串行（$2^k-1$ 次）；
- 通过 `alg.trunc` 可读取 `algmult.trunc`。

### 构建 IF 的入口

```julia
hybriddynamics(lattice::AbstractADTLattice, corr, bs::AdditiveHyb; kwargs...)
hybriddynamics(lattice::AbstractADTLattice, corr, bs::AdditiveHyb, alg::InfluenceFunctionalAlgorithm)
hybriddynamics(lattice::AbstractPTLattice, corr, hyb::GeneralHybStyle, alg::TranslationInvariantIF)
hybriddynamics!(gmps, lattice, corr, hyb[, alg]; kwargs...)

hybriddynamics_naive(lattice, corr, bs::AdditiveHyb; trunc=DefaultITruncation)
hybriddynamics_naive!(gmps, lattice, corr, bs::AdditiveHyb; trunc=DefaultITruncation)
```

- `hybriddynamics(lattice, corr, hyb, trunc=trunc)`：对角耦合 + 部分 IF（默认 `PartialIF(trunc=DefaultITruncation)`），返回 ADT；
- `hybriddynamics(lattice, corr, hyb, alg)`：按 `alg` 构建（对角也可用 `TranslationInvariantIF`），非对角必须用 `TranslationInvariantIF`，返回 ProcessTensor；
- `hybriddynamics_naive`：$O(N^2)$ 门操作的朴素参考实现（仅对角，用于验证）。

### 单行部分 IF

```julia
partialif(lattice::AbstractADTLattice, rowind::ContourIndex, corr, hyb::AdditiveHyb)
partialif_naive(lattice::AbstractADTLattice, rowind::ContourIndex, corr, hyb::AdditiveHyb; trunc=DefaultITruncation)
```

构造 IF 中固定"行"（时间点 `rowind`）的单部分 MPS。`partialif` 精确版（键维 2），`partialif_naive` 直接指数化并逐步应用。

### 影响算符（底层工具）

```julia
influenceoperator(lattice, corr, hyb::GeneralHybStyle; algexpan=OverDeterminedProny())
influenceoperatorexponential(lattice, corr, dt::Real, hyb::GeneralHybStyle, alg; algexpan=OverDeterminedProny())
differentialinfluencefunctional(lattice, corr, dt::Real, hyb::GeneralHybStyle, alg, algmult::DMRGAlgorithm; algexpan=OverDeterminedProny())
```

- `influenceoperator`：返回有效哈密顿量 $H_{\text{eff}}$ 的分支 MPO（`mpo₁,mpo₂,mpo₃,mpo₄`）；
- `influenceoperatorexponential`：返回时间步进后的 MPO（`FirstOrderStepper` 返回 4 个，`ComplexStepper` 返回 8 个）；
- `differentialinfluencefunctional`：把步进 MPO 相乘得到最小时间步的 IF（`DMRGAlgorithm` 指定压缩）。

---

## 12. 边界条件与初态

```julia
boundarycondition(x::ADT, lattice; kwargs...)          # 非就地（copy）
boundarycondition!(x::ADT, lattice; ρ₀::VecOrMat=ones(lattice.d))
initialstate!(x::ADT, lattice::RealADTLattice1Order, ρ0::AbstractVecOrMat)
initialstate!(x::ProcessTensor, lattice::RealPTLattice1Order, ρ0::AbstractMatrix)
```

- 实时间 ADT：`boundarycondition!` 施加单位边界算符并调用 `initialstate!` 设置初态密度矩阵 `ρ₀`（默认最大混合 $\mathbb{1}/d$，也可传对角向量）；
- 实时间 PT：`initialstate!` 把初态 $\rho_0$ 吸收进前两个位点（带 SVD 压缩，`DefaultIntegrationTruncation`）；
- 虚时间 ADT：`boundarycondition!` 施加单位算符（对应周期/开边界配分函数）；
- 混合轮廓：三处轮廓连接点施加单位算符。

---

## 13. 系统模型 Models

### ImpurityHamiltonian

```julia
struct ImpurityHamiltonian{M<:AbstractMatrix} <: AbstractImpurityOperator
    m::M
end
ImpurityHamiltonian(m::AbstractMatrix)
ImpurityHamiltonian(d::Int)     # 零矩阵 d×d
```

幺正系统：每个时间步施加传播子 $e^{\mp i H\delta t}$（虚时间 $e^{-H\delta\tau}$）。

### ImpurityLindbladian

```julia
struct ImpurityLindbladian <: AbstractImpurityOperator
    m::Array{ComplexF64, 4}
end
ImpurityLindbladian(L::LindbladOperator)
ImpurityLindbladian(H::AbstractMatrix, jumpops::Vector{<:AbstractMatrix})
ImpurityLindbladian(h::ImpurityHamiltonian)
ImpurityLindbladian(d::Int)
```

含耗散系统：每个时间步施加 Lindblad 传播子 $e^{\mathcal{L}\delta t}$（`LindbladOperator` 由 `ImpurityModelBase` 提供）。

### sysdynamics / sysdynamics!

```julia
sysdynamics(lattice, model::AbstractImpurityOperator; trunc::TruncationScheme=DefaultKTruncation)
sysdynamics(lattice, model, op::ContourOperator; trunc=DefaultKTruncation)
sysdynamics!(mps, lattice, model[, op]; trunc=DefaultKTruncation)
```

构建裸系统动力学张量 `mpsK`（不含影响泛函）：
- `op` 关键字：在系统传播子中插入算符（用于关联函数）；
- 每个站点应用传播子后做 `canonicalize!`（`trunc` 控制截断，默认 `DefaultKTruncation` D=1000）。

---

## 14. 观测量 Observables

### environments

```julia
environments(lattice::AbstractADTLattice, A::ADT, B::Vararg{ADT})
environments(lattice::RealPTLattice, A::ProcessTensor, B::Vararg{ProcessTensor}; ρ₀::AbstractMatrix=_eye(phydim(lattice)))
environments(lattice::ImagPTLattice, A::ProcessTensor, B::Vararg{ProcessTensor})
environments(lattice::MixedPTLattice, A::ProcessTensor, B::Vararg{ProcessTensor})
```

预计算张量网络的全环境（转移矩阵缩并缓存），返回缓存对象（`ADTExpectationCache` / `PTExpectationCache`）。实时间 PT 需要初态密度矩阵 `ρ₀`。常见调用：
- ADT 实时间：`environments(lattice, mpsK, mpsI)`；
- PT 实时间：`environments(lattice, mps, ρ₀=ρimp)`；
- 虚时间：`environments(lattice, mpsK, mpsI)` 或 `environments(lattice, mps)`。

### expectationvalue / expectation

```julia
expectationvalue(m::ADTTerm, cache::ADTExpectationCache)
expectationvalue(m::ContourOperator, cache::PTExpectationCache)
expectationvalue(m::AbstractFockTerm, cache::PTExpectationCache)

expectation(m::ADTTerm, cache::ADTExpectationCache)
expectation(m::AbstractFockTerm, cache::PTExpectationCache)
```

- `expectationvalue`：归一化期望值 $\langle \hat{O}\rangle = \mathrm{expectation}/Z$；
- `expectation`：未归一化的缩并值；
- `ADTTerm` 用于 ADT 缓存，`ContourOperator` / `FockTerm` 用于 PT 缓存。

### Zvalue / Zvalue2

```julia
Zvalue(cache)    # 配分函数 Z（缓存左端与右端边界缩并）
Zvalue2(cache)   # 另一种边界组合（PT 缓存）
```

### TransferMatrix

```julia
TransferMatrix(states::Vararg{ADT})                  # → ADTTransferMatrix
TransferMatrix(j::Int, states::Vararg{ADT})          # 单点转移矩阵
TransferMatrix(lattice::AbstractPTLattice, states::Vararg{ProcessTensor}; scaling=scaling(states...))
```

转移矩阵对象（`ADTTransferMatrix` / `PTTransferMatrix`），支持 `left * m`、`m * right` 形式的多重转移缩并，是 `environments` 的基础。PT 版本对实时间格点自动把 2 个位点合并为一步。

---

## 附：默认值速查（`src/defaults.jl`）

| 常量 | 值 | 用途 |
|---|---|---|
| `DefaultTruncation` | `truncdimcutoff(D=100, ϵ=1e-14)` | 通用默认 |
| `DefaultITruncation` | `truncdimcutoff(D=200, ϵ=1e-10)` | 构建 IF / mult 默认 |
| `DefaultKTruncation` | `truncdimcutoff(D=1000, ϵ=1e-10)` | 系统动力学默认 |
| `DefaultIntegrationTruncation` | `truncdimcutoff(D=10000, ϵ=1e-12)` | 初态吸收 |
| `DefaultMPOTruncation` | `truncdimcutoff(D=10000, ϵ=1e-12)` | MPO 压缩 |
| `DefaultMultAlg` | `DMRGMult1(DefaultITruncation)` | `mult` 默认算法 |

---

*配合 `docs/documentation.md`（概念与示例）阅读。*
