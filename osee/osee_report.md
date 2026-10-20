# 玻色杂质影响泛函的 OSEE：解析计算与 TEMPO 数值验证

## 1. 问题与结论

**问题**：若杂质也是玻色子、bath 为高斯玻色 bath，Feynman–Vernon 影响泛函（IF）在指数上为二次型，整体等价于高斯算符。能否由 hybridization function Δ 快速求出 IF 的 operator-space entanglement entropy（OSEE），并与 TEMPO 构造的 PT/ADT 沿时间 bipartition 的纠缠熵对比？

**结论**：

1. TEMPO 的 IF 为 `F = exp(+H)`，其中 `H = Σᵢⱼ Δᵢⱼ a†ᵢaⱼ`（Δ 来自 `correlationfunction(bath, lattice)`）。
2. 用 ImpurityModelBase 构造同一有限 Fock 截断 `d` 下的 quadratic 模型，直接对 `F = exp(H)` 做算符 SVD，得到的 OSEE 与 TEMPO PT/ADT 的 bond 熵一致：
   - PT：d=3 误差 ~1e-6，d=6 误差 ~1e-5；
   - ADT（对角耦合）：k=5 误差 ~0.9%，k=7 误差 ~0.2%（残差来自 XTRGIF 截断）。
3. 弱耦合下可由 Δ 的非对角块快速估算：
   `S_l ≈ 2 Σₖ sₖ²(1 − ln sₖ²)`，`sₖ = svdvals(Δ[1:l, l+1:N])`，误差约几个百分点。
4. **不存在仅由 Δ 决定、与 Fock 截断 d 无关的 OSEE**。固定 δτ 时 d→∞ 的极限有限（由 Δ 的最小正本征值 λ_min 决定，量级 ~ln(1/λ_min)），但收敛尺度 d ~ 1/λ_min；δτ→0 时该极限对数发散。实际有物理意义且与 TEMPO 计算一致的是**固定 d** 的 OSEE。

---

## 2. 模型与参数

Ohmic 谱：

```
J(ω) = 2α ω exp(−ω/ωc),   α=0.1, ωc=2.0
```

默认测试参数（若不特别说明）：

```
N = 4, δτ = 0.1, β = N δτ = 0.4, d = 3 (Fock 截断)
```

对应的 Δ 矩阵（`correlationfunction(bath, PTLattice(...)).data`，虚时）：

```
 0.006369     0.00579114   0.00628409   0.00684059
 0.00684059   0.006369     0.00579114   0.00628409
 0.00628409   0.00684059   0.006369     0.00579114
 0.00579114   0.00628409   0.00684059   0.006369
```

本征值：

```
[2.1367e-5, 8.49135e-5 ± 0.0010494456im, 0.0252848224]
```

Δ 是循环矩阵但**非对称**（Δ₁₂ ≠ Δ₂₁）。最小正本征值 `λ_min ≈ 2.1e-5` 决定了 d 收敛尺度。

---

## 3. OSEE 的定义与计算

把算符 `F` 向量化为 doubled Hilbert space 中的（未归一化）态
`|F⟩ = (F ⊗ I)|I⟩`。对时间 cut `l`（前 l 个时间点 vs 其余），把系数张量
`F[n₁…n_N, m₁…m_N]` 按 `{n_{≤l}, m_{≤l}}` 与 `{n_{>l}, m_{>l}}` 重排，
做 SVD 得奇异值 `s_k`，令 `p_k = s_k²/Σs²`，则

```
S_l = −Σ_k p_k ln p_k
```

Julia 实现（与 ImpurityModelBase 的 MSB 站点顺序一致）：

```julia
function op_entropy(F, d, N, l)
    R = reshape(F, ntuple(_ -> d, 2N)...)
    rowdims = [N - s + 1 for s in 1:N]
    coldims = [N + (N - s + 1) for s in 1:N]
    perm = vcat(vcat(rowdims[1:l], coldims[1:l]),
                vcat(rowdims[l+1:N], coldims[l+1:N]))
    R2 = permutedims(R, perm)
    Mmat = reshape(R2, d^(2l), d^(2(N - l)))
    sv = svdvals(Mmat)
    p = sv .^ 2; p ./= sum(p)
    p = p[p .> 1e-18]
    return -sum(p .* log.(p))
end
```

精确 quadratic 模型（PT 耦合 `a†`）：

```julia
terms = [adaga(i, j, coeff=Δ[i, j]) for i in 1:N for j in 1:N if Δ[i, j] != 0]
H = bosonoperator(quadratichamiltonian(N, terms); d=d)
F = exp(H)                       # 注意符号为 +H
S = [op_entropy(F, d, N, l) for l in 1:N-1]
```

ADT 对角耦合（`AdditiveHyb([0,1,…,d-1])` = n̂）对应
`F = exp(Σᵢⱼ Δᵢⱼ nᵢnⱼ)`，可用 dense 的 `kron` 直接构造。

---

## 4. 数值验证

### 4.1 PT 路径（耦合 NonDiagonalHyb(bosonadagoperator(d))）

| d | 精确 exp(+H) OSEE (l=1,2,3) | TEMPO PT bond 熵 | 误差 |
|---|---|---|---|
| 3 | 0.0024671934, 0.0031871446, 0.0024671934 | 0.0024663474, 0.0031874827, 0.0024663474 | ~1e-6 |
| 6 | 0.0131641397, 0.0168897048, 0.0131641397 | 0.0131500785, 0.0168785472, 0.0131500785 | ~1e-5 |

PT 设置：`XTRGIF(k=5, fast=true, D=200, ϵ=1e-14)`，`OverDeterminedProny(n=3, tol=1e-12)`；
bond dimension d=3 时 59，d=6 时 146。

**符号确认**：`exp(+H)` 与 PT 的误差 0.0083，`exp(−H)` 为 0.070；且 d=6 的熵匹配到 1e-5，确认 TEMPO 的 IF 是 `exp(+H)`。

### 4.2 精确 OSEE 随 d 的变化（PT 耦合）

| d | S(l=1,2,3) | ‖F‖ |
|---|---|---|
| 3 | 0.0024671934, 0.0031871446, 0.0024671934 | 9.237728 |
| 4 | 0.0051777767, 0.0066711007, 0.0051777767 | 16.645099 |
| 5 | 0.0087567291, 0.0112572609, 0.0087567291 | 26.368013 |
| 6 | 0.0131641397, 0.0168897048, 0.0131641397 | 38.507255 |
| 8 | 0.02436699, 0.0311515263, 0.02436699 | 70.475616 |

小 d 区间近似 ∝ d²。

### 4.3 ADT 路径（对角耦合 n̂）

`AdditiveHyb(collect(0.0:d-1))`，`ImagADTLattice(N=4, δτ=0.1, d=3)`：

| 模型 | l=1 | l=2 | l=3 |
|---|---|---|---|
| 精确 exp(+ΣΔᵢⱼnᵢnⱼ) | 0.0020616 | 0.0026610 | 0.0020616 |
| 精确 exp(−ΣΔᵢⱼnᵢnⱼ) | 0.0019160 | 0.0024735 | 0.0019160 |
| ADT k=5 | 0.0020430 | 0.0026382 | 0.0020430 |
| ADT k=7 | 0.0020569 | 0.0026561 | 0.0020569 |

k=5→7 误差由 0.9% 降到 0.2%，说明残差来自 `(I + H/2ᵏ)^{2ᵏ}` 的 XTRGIF 截断，而非公式错误。ADT 的 `.s` 首末 bond 为真空边界（S=0）。

### 4.4 Δ 对角项的作用

| 模型（PT，d=3） | l=1 | l=2 | l=3 |
|---|---|---|---|
| exp(Σ**全部** Δᵢⱼa†ᵢaⱼ) | 0.00246719 | 0.00318714 | 0.00246719 |
| exp(Σ**i≠j** Δᵢⱼa†ᵢaⱼ) | 0.00245783 | 0.00317509 | 0.00245783 |
| TEMPO PT | 0.00246635 | 0.00318748 | 0.00246635 |

**含对角项才是正确对照**（TEMPO 与“全部项”差 1e-6，与“仅非对角”差约 0.4%）。

物理解释：Δᵢᵢ 只作用在同一时间点，`ΣᵢΔᵢᵢa†ᵢaᵢ` 在时间 cut 下是 `O_L ⊗ O_R` 的局域非幺正滤波，不产生跨时间纠缠，但会改变 Schmidt 谱，因此仍改变熵的数值（弱耦合下约 0.4%）。所以：

- 问“纠缠从哪来”：只有非对角 Δᵢⱼ（i≠j）；
- 问“熵的数值”：必须保留对角项；只在快速近似公式中可以丢掉。

---

## 5. 由 Δ 快速估算 OSEE

### 5.1 弱耦合 SVD 公式

```
S_l ≈ 2 Σₖ sₖ² (1 − ln sₖ²),   sₖ = svdvals(Δ[1:l, l+1:N])
```

结果（N=4，Δ 同上）：

| l | 奇异值 | S_fast | d=3 精确 |
|---|---|---|---|
| 1 | 0.01094627 | 0.0024034874 | 0.0024671934 |
| 2 | 0.01262182, 9.9e-6 | 0.0031048554 | 0.0031871446 |
| 3 | 0.01094627 | 0.0024034874 | 0.0024671934 |

误差约 2.6%。来源：非对角项贡献 `~|Δ_off|² ln(1/|Δ_off|²)`，对角项只是高阶局域修正。

### 5.2 高斯协方差公式（d→∞ 极限）

归一化 IF 可视为无限维高斯算符。取 `M = exp(−Δ)`（`exp(+Δ)` 只差每个模的粒子-空穴反射，时间 cut 熵相同），构造

```
W = [[0, M], [Mᵀ, 0]]
G = (I − W W)^{-1}          # ⟨c c†⟩
K = G W                     # ⟨c c⟩
σ_xx = G + K − ½I
σ_pp = G − K − ½I
Ω = [[0, I], [−I, 0]]
ν = positive eigvals(i Ω σ)
S = Σ_ν [(ν+½)ln(ν+½) − (ν−½)ln(ν−½)]
```

注意：**必须减 ½I**，且只取正辛本征值（早期漏掉这两点导致单对自检结果翻倍）。

单对压缩态 `|ψ⟩ = Σ sⁿ|n,n⟩` 自检：

| s | 公式 | 解析 `−ln(1−s²) − s²/(1−s²)ln s²` |
|---|---|---|
| 0.2 | 0.1749418206 | 0.1749418206 |
| 0.5 | 0.7497801928 | 0.7497801928 |
| 0.8 | 1.8150505411 | 1.8150505411 |

真空整体系熵为 0。

对本模型（N=4，M=exp(−Δ)）：

```
l=1  S_∞ = 6.0336782129
l=2  S_∞ = 6.0954509270
l=3  S_∞ = 6.0336786073
```

这**不是**与有限 d 无关的物理值，而是 d→∞ 的极限。`exp(+Δ)` 直接用该公式会因 `I−W²` 不定（奇异值 >1）数值失败。

---

## 6. d→∞ 是否发散？

固定 δτ 时**不发散**，而是收敛到上面高斯公式给出的有限值。用 N=2、按总粒子数分块精确计算 `exp(H)`（H 守恒总粒子数，每块 ≤ d 维），Δ 的本征值为 `[2.14e-5, 0.02528]`：

| d | 4 | 8 | 16 | 32 | 64 | 128 | 192 | 256 |
|---|---|---|---|---|---|---|---|---|
| OSEE | 2.7703 | 4.1469 | 5.4877 | 6.6590 | 7.1840 | 6.8477 | 6.6046 | 6.4497 |

峰值在 d≈64，之后下降并趋向高斯极限 **6.304**。所以有限 d 的小熵（0.002）与高斯极限（6.03）是同一条曲线的两端，不是两个无关的数。

需要区分的情况：

1. **δτ 固定、λ_min>0**：d→∞ 极限有限，但收敛尺度 `d ~ 1/λ_min`；本模型 λ_min≈2.1e-5，故 d 要到 ~5×10⁴ 才饱和。在 d 为几十到几千时表现为增长（小 d 近似 ∝d²），看起来像发散。
2. **δτ→0**：λ_min→0，极限值 `~ln(1/λ_min)` **对数发散**，OSEE 不是连续极限下的普适物理量。
3. **Δ 有严格零本征值**：对应模式完全混合，`S ~ ln d`，d→∞ 真发散。

---

## 7. 复现脚本要点

环境：

```bash
JULIA_LIBSTDCXX_PROBE=0 julia --startup-file=no <script>.jl
```

核心流程：

```julia
using TEMPO, ImpurityModelBase, LinearAlgebra

omic_spectrum(w, α, wc) = 2α * w * exp(-w / wc)
N, α, wc = 4, 0.1, 2.0
δτ = 0.1; β = N * δτ; d = 3

spec  = spectrum(w -> omic_spectrum(w, α, wc), lb=0, ub=wc)
bath  = bosonicbath(spec, β=β)
Δ     = Matrix(correlationfunction(bath, PTLattice(N=N, δτ=δτ, d=d, contour=:imag)).data)

# 精确 quadratic OSEE（PT 耦合）
terms = [adaga(i, j, coeff=Δ[i, j]) for i in 1:N for j in 1:N if Δ[i, j] != 0]
H = bosonoperator(quadratichamiltonian(N, terms); d=d)
S_exact = [op_entropy(exp(H), d, N, l) for l in 1:N-1]

# TEMPO PT
trunc  = truncdimcutoff(D=200, ϵ=1e-14, add_back=0)
alg    = XTRGIF(k=5, fast=true,
                alg_mult=TEMPO.DMRG1(trunc=trunc),
                alg_expan=OverDeterminedProny(n=3, tol=1e-12))
pt     = hybriddynamics(PTLattice(N=N, δτ=δτ, d=d, contour=:imag),
                        correlationfunction(bath, PTLattice(N=N, δτ=δτ, d=d, contour=:imag)),
                        NonDiagonalHyb(bosonadagoperator(d=d)), alg)
S_pt   = 每个 bond 的归一化奇异值谱熵
```

ADT 路径只需把 lattice 换成 `ADTLattice(N=N, δτ=δτ, d=d, contour=:imag)`，耦合换成
`AdditiveHyb(collect(0.0:d-1))`。

同目录下的测试脚本：

- `osee_fast.jl`：精确 quadratic OSEE（d 扫描）、TEMPO PT 对比、弱耦合 SVD 快速公式；
- `adt_check.jl`：ADT 对角耦合路径与精确 `exp(ΣΔᵢⱼnᵢnⱼ)` 对比（k=5/7）；
- `gauss_check.jl`：高斯协方差公式的单对自检与 `M=exp(±Δ)` 极限值；
- `osee_dlimit.jl`：N=2 按总粒子数分块扫描 d，验证 OSEE 收敛到高斯极限而非发散。

运行方式：

```bash
JULIA_LIBSTDCXX_PROBE=0 julia --startup-file=no osee_dlimit.jl
```

---

## 8. 实践建议

1. **要精确、要与 TEMPO 对应**：固定 Fock 截断 d，构造 `exp(+ΣΔᵢⱼa†ᵢaⱼ)` 做算符 SVD。这是与 PT/ADT bond 熵严格一致的量。
2. **要快速估算**：用非对角块 SVD 公式，弱耦合误差约几个百分点；不要用它代替精确值。
3. **不要**试图用 `M=exp(+Δ)` 直接套高斯协方差公式；如需 d→∞ 极限，用 `M=exp(−Δ)`，并记住该值由 λ_min 决定、对 δτ 敏感，不是普适物理量。
4. 对角项 Δᵢᵢ 不产生时间纠缠，但会改变熵的数值；精确计算中必须保留。
