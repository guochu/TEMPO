# API 参考

本页按功能分类列出 `TEMPO` 的导出符号文档（由 docstring 自动生成）。使用方式：

```julia
using TEMPO, ImpurityModelBase, LinearAlgebra
```

说明：
- `branch`、`index` 由 `QuAPI` 提供，`scalartype`、`space_l`、`space_r` 由 `TensorOperations` 提供，经 TEMPO 重导出；
- 指数展开函数（`OverDeterminedProny`、`DeterminedProny`、`exponential_expansion`、`expansion_error`）由 `ExpExp` 包提供并重导出，见文末[指数展开（`ExpExp`，重导出）](@ref)小节的签名速查；
- `phydim`、`spectrum`、`bosonicbath` 等在 `ImpurityModelBase` 中定义。

## 截断方案

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["auxiliary/auxiliary.jl", "auxiliary/truncation.jl", "auxiliary/distance.jl",
         "conversions.jl", "defaults.jl"]
```

## 张量分解与矩阵操作

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["auxiliary/matrixalgebra.jl", "auxiliary/tensorops.jl", "auxiliary/distance.jl"]
```

## MPS/DMRG 乘法算法

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["auxiliary/mpsalgs.jl"]
```

## 轮廓索引

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["contourindices.jl"]
```

## MPO 哈密顿量

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["mpohamiltonian/def.jl", "mpohamiltonian/abstractmpotensor.jl",
         "mpohamiltonian/sparsempotensor.jl", "mpohamiltonian/schurmpotensor.jl",
         "mpohamiltonian/mpohamiltonian.jl",
         "mpohamiltonian/schurmpo/schurmpo.jl", "mpohamiltonian/schurmpo/longrange.jl",
         "mpohamiltonian/schurmpo/exponentialdecay.jl", "mpohamiltonian/schurmpo/generaldecay.jl",
         "mpohamiltonian/schurmpo/w1w2.jl"]
```

## 增强密度张量 ADT

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["adt/def.jl", "adt/adt.jl", "adt/abstractdefs.jl", "adt/orth.jl",
         "adt/mult/mult.jl", "adt/mult/svdmult.jl", "adt/mult/iterativemult.jl",
         "adt/integrate.jl", "adt/util.jl", "adt/linalg.jl"]
```

## 过程张量 ProcessTensor

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["pt/def.jl", "pt/orth.jl", "pt/mult/mult.jl", "pt/mult/svdmult.jl",
         "pt/mult/iterativemult.jl", "pt/util.jl", "pt/linalg.jl"]
```

## 局域算符项（ADTTerm / FockTerm）

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["adtterms.jl", "fockterms.jl", "contouroperators.jl"]
```

## 格点 Lattice

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["adtlattices/adtlattices.jl", "adtlattices/fockordering.jl",
         "adtlattices/realtime.jl", "adtlattices/imaginarytime.jl", "adtlattices/mixedtime.jl",
         "ptlattices/ptlattices.jl", "ptlattices/realtime.jl", "ptlattices/imaginarytime.jl",
         "ptlattices/mixedtime.jl", "ptlattices/integrate.jl"]
```

## 浴关联函数

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["correlationfunction.jl"]
```

## 影响泛函 IF

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["influencefunctional/influencefunctional.jl",
         "influencefunctional/partialif/partialif.jl", "influencefunctional/partialif/util.jl",
         "influencefunctional/partialif/realtime.jl", "influencefunctional/partialif/imaginarytime.jl",
         "influencefunctional/partialif/mixedtime.jl",
         "influencefunctional/ptpartialif/ptpartialif.jl", "influencefunctional/ptpartialif/util.jl",
         "influencefunctional/ptpartialif/realtime.jl", "influencefunctional/ptpartialif/imaginarytime.jl",
         "influencefunctional/ptpartialif/mixedtime.jl",
         "influencefunctional/ttiif/ttiif.jl",
         "influencefunctional/ttiif/adt/adt.jl", "influencefunctional/ttiif/adt/imag.jl",
         "influencefunctional/ttiif/adt/real.jl",
         "influencefunctional/ttiif/pt/pt.jl", "influencefunctional/ttiif/pt/imag.jl",
         "influencefunctional/ttiif/pt/real.jl"]
```

## 含时影响泛函与 BEC 影响泛函

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["tdinfluencefunctional/tdinfluencefunctional.jl",
         "tdinfluencefunctional/partialif/partialif.jl", "tdinfluencefunctional/ttiif/ttiif.jl",
         "becinfluencefunctional/becinfluencefunctional.jl"]
```

## 边界条件与初态

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["boundarycondition.jl"]
```

## 系统模型 Models

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["models/def.jl", "models/models.jl", "models/dissipative.jl",
         "models/unitary/unitary.jl", "models/unitary/adt.jl", "models/unitary/pt.jl"]
```

## 观测量 Observables

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["observables/adt/adt.jl", "observables/adt/envs.jl", "observables/adt/transfer.jl",
         "observables/adt/gf.jl",
         "observables/pt/pt.jl", "observables/pt/envs.jl", "observables/pt/transfer.jl",
         "observables/pt/mixedtransfer.jl",
         "observables/observables.jl",
         "observables/correlations.jl", "observables/heatcurrents.jl"]
```

## 指数展开（`ExpExp`，重导出）

指数展开将数据 $f(x)$（$x = 1,\dots,N$）拟合为指数和 $f(x) \approx \sum_i \alpha_i \beta_i^{x}$，用于把混合化函数表示为 MPO 可演化的形式。

```julia
OverDeterminedProny(; n::Int=10, tol::Real=1.0e-8, verbosity::Int=1, stepsize=1)
DeterminedProny(; n::Int=10, tol::Real=1.0e-8, verbosity::Int=1, stepsize=1)
#   n         最大展开项数（迭代上限）
#   tol       相对误差阈值（tol*norm(f)），达到即提前收敛
#   verbosity ≥1 打印不收敛警告，≥2 打印收敛信息
#   stepsize  均匀采样步长；nothing 时自动选择

exponential_expansion(f::Vector{<:Number}, [L::Int,] alg=OverDeterminedProny())
#   → (αs, βs)：满足 f(x) ≈ Σᵢ αᵢ βᵢˣ；传入函数时自动采样 k=1:L

expansion_error(f, p) / expansion_error(f, coeffs, alphas)
#   展开在采样点上的 2-范数误差
```

## 默认值速查（`src/defaults.jl`）

| 常量 | 值 | 用途 |
|---|---|---|
| `DefaultTruncation` | `truncdimcutoff(D=100, ϵ=1e-14)` | 通用默认 |
| `DefaultITruncation` | `truncdimcutoff(D=200, ϵ=1e-10)` | 构建 IF / mult 默认 |
| `DefaultKTruncation` | `truncdimcutoff(D=1000, ϵ=1e-10)` | 系统动力学默认 |
| `DefaultIntegrationTruncation` | `truncdimcutoff(D=10000, ϵ=1e-12)` | 初态吸收 |
| `DefaultMPOTruncation` | `truncdimcutoff(D=10000, ϵ=1e-12)` | MPO 压缩 |
| `DefaultMultAlg` | `DMRGMult1(DefaultITruncation)` | `mult` 默认算法 |
