# 快速上手

本页给出四类典型问题的完整计算流程，代码可直接照抄运行。选择哪一段取决于**耦合形式**与**想要的观测量**：

| 问题 | 耦合 | 格点 | 示例 |
|---|---|---|---|
| [标准自旋玻色子模型（对角耦合，Rabi 型，ADT）](@ref) | 对角 | `ADTLattice` | `docs/tutorials/spinboson/rabitype.jl` |
| [JC 型非对角耦合（PT 框架）](@ref) | 非对角（共轭对） | `PTLattice` | `docs/tutorials/spinboson/jctype.jl` |
| [玻色杂质：虚时间演化与 Matsubara 关联函数](@ref) | 对角或非对角 | `ADTLattice`/`PTLattice`，`contour=:imag` | `benchmark/independentbosons.jl` |
| [含时系统-浴耦合](@ref) | 任意 | 任意 | `src/tdinfluencefunctional/` |

所有计算都遵循同一模板：

```text
1. trunc  = 截断方案
2. lattice = ADTLattice/PTLattice（轮廓 + 步数 + 步长）
3. hyb    = AdditiveHyb / NonDiagonalHyb / ...（耦合算符）
4. spec/bath = spectrum + bosonicbath（浴谱与温度）
5. corr   = correlationfunction(bath, lattice)
6. mpsI   = hybriddynamics(lattice, corr, hyb[, alg])    # 影响泛函
7. mpsK   = sysdynamics(lattice, model) + boundarycondition!（系统动力学 + 初态）
8. cache  = environments(...)；expectationvalue(...)      # 观测量
```

## 标准自旋玻色子模型（对角耦合，Rabi 型，ADT）

系统是二能级自旋，$\hat{A}=\hat{\sigma}_x/2$（对角耦合），浴谱取 sub-Ohmic 谱。

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
hyb  = AdditiveHyb(zdiag)                       # 对角耦合（注意传对角元向量）
spec = spectrum(w -> 2π*α * w^s * wc^(1-s), lb=0, ub=wc)   # sub-Ohmic 谱
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

# 4. 构建影响泛函（IF），得到一个 ADT（MPS）
mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)

# 5. 构建裸系统动力学并施加边界条件（初态）
model = ImpurityHamiltonian(Δ .* x)             # 系统哈密顿量 H_S
mpsK  = sysdynamics(lattice, model, trunc=trunc)
mpsK  = boundarycondition!(mpsK, lattice, ρ₀=ρimp)

# 6. 预计算环境，用于观测量
cache = environments(lattice, mpsK, mpsI)

# 7. 计算局部观测量（如 ⟨σ_z(t)⟩）
obs = ComplexF64[]
for i in 1:Nt
    pos = index(lattice, i, branch=:+)
    m   = ADTTerm(pos, zdiag)
    v   = expectationvalue(m, cache)            # 已除以 Z，得到归一化期望值
    push!(obs, v)
end
```

!!! note "对角耦合也可以用 XTRG-IF"
    上面的 `hybriddynamics(lattice, corr, hyb, trunc=trunc)` 走 `PartialIF` 路径。对更大的局域维数（如 spin-1，d=3）或需要更精确的压缩时，用 `XTRGIF` 算法（见[使用手册](@ref)）通常更快、内存更省，见[实践指南](@ref)的讨论。

!!! tip "ADT 路径也能测非对角算符：算符插入"
    除了 `ADTTerm`，ADT 格点上还可以通过**向系统动力学中插入算符**来测非对角算符与两点关联函数（虚时间测试见 `test/models/rabimodel.jl`）：

    ```julia
    # 单点算符（可以是任意矩阵，含非对角元）
    ct = ContourOperator(ContourIndex(1), op)          # op 任意矩阵
    mpsK = sysdynamics(lattice, model, ct, trunc=trunc)

    # 两点关联 ⟨op2(t_i) op1(t_1)⟩：两个位置各插一个算符
    ct = ContourOperator([ContourIndex(i), ContourIndex(1)], [op2, op1])
    mpsK = sysdynamics(lattice, model, ct, trunc=trunc)

    mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)  # 实时间需要 ρ₀
    mps2 = mult!(mpsK, mpsI, trunc=trunc)
    v = integrate(mps2) / integrate(mps)               # 已归一化
    ```

    这种方式每次插入都需重建 `mpsK` 并重新乘 IF，适合少量算符/关联函数；批量测单点对角量时 `ADTTerm` + `environments` 缓存更高效。ADT 上同样可以用 `ADTTerm((pos2, pos1), (v2, v1))` 的多点形式测对角两点关联（`apply!` + `integrate(mps2)/integrate(mps)`）。

## JC 型非对角耦合（PT 框架）

系统通过 $\hat{A}=\hat{\sigma}_-/2$（共轭对 $\hat{A}^\dagger,\hat{A}$）耦合到浴，这是文献中的**非对角耦合**情形，必须使用 PT 框架与平移不变影响泛函算法。

```julia
lattice = PTLattice(N=Nt, δt=δt, contour=:real)   # 注意是 PTLattice

hyb  = NonDiagonalHyb(sp)                        # 非对角耦合：op*a + op'*a'
spec = spectrum(w -> 2π*α * w^s * wc^(1-s), lb=0, ub=wc)
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

# IF 构建算法：平移不变 IF（XTRG 式），
# 使用 DMRG 型 MPO-MPO 乘法 + Prony 指数展开
algmult  = DMRGMult1(trunc, maxiter=10)
algexpan = OverDeterminedProny(n=n, tol=1.0e-8, verbosity=2)
alg = XTRGIF(k=k, fast=true, algmult=algmult, algexpan=algexpan, verbosity=2)

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

!!! warning "测非对角观量的两条路线"
    ADT 路径的 `ADTTerm` 只支持对角元向量，要测 $\langle\sigma_x\rangle$ 等非对角量，可以：① 改用 `PTLattice` + `NonDiagonalHyb`（对实对角算符与 `AdditiveHyb` 物理等价），用 `ContourOperator` 批量测（推荐，见上一节示例）；② 保持在 ADT 路径用**算符插入**（见上文 tip）。详见[实践指南](@ref)。

## 玻色杂质：虚时间演化与 Matsubara 关联函数

此时杂质本身是玻色模（局域 Hilbert 空间截断为 `d`），在虚时间轮廓上计算 Matsubara Green 函数 $\langle \mathcal{T}_\tau \hat{a}(\tau)\hat{a}^\dagger(0)\rangle$。

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

对**非对角耦合**的玻色杂质（如 `benchmark/bosonicimpurity.jl`）：

```julia
lattice = PTLattice(N=N, δt=δt, d=d, contour=:real)
hyb  = NonDiagonalHyb(a')
alg  = XTRGIF(k=5, fast=true,
                             algmult=DMRGMult1(trunc, initguess=:rand),
                             algexpan=OverDeterminedProny(n=20, tol=1.0e-8))
mpsI = hybriddynamics(lattice, corr, hyb, alg)
```

## 含时系统-浴耦合

对含时耦合（见 `src/tdinfluencefunctional/`），使用含时版本的混合化样式：`AdditiveTdHyb`、`NonAdditiveTdHyb`、`NonDiagonalTdHyb`。它们接受一个函数 `f(t)` 描述耦合强度的时间依赖：

```julia
hyb = NonDiagonalTdHyb(op, t -> f(t))   # 或 AdditiveTdHyb / NonAdditiveTdHyb
```
