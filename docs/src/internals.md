# 实现细节

本文档深入介绍 `TEMPO` 工具包的内部实现：数据结构、算法流程、张量约定与数值细节，与文献（Phys. Rev. B 114, 125413 (2026)）中的公式一一对应。接口与用法见[快速上手](@ref)与[使用手册](@ref)。

## 1. 总览与模块结构

工具包在 `ImpurityModelBase`（系统模型定义）与 `QuAPI`（基础张量网络类型 `Dense1DTN`、`DenseMPSTensor`、`ContourIndex`、`branch`/`index` 等）之上实现 TEMPO 算法。加载顺序见 `src/TEMPO.jl`：

```
auxiliary ──→ defaults ──→ mpohamiltonian ──→ contourindices
      └─→ adt ──→ pt ──→ conversions ──→ adtterms / fockterms
      └─→ adtlattices ──→ ptlattices ──→ contouroperators
      └─→ correlationfunction ──→ influencefunctional ──→ tdinfluencefunctional
      └─→ boundarycondition ──→ models ──→ observables
```

各模块职责：

| 模块 | 职责 |
|---|---|
| `auxiliary/` | 截断方案、指数展开、MPS 算法（`DMRGMult1`、`SVDCompression`）、`transfermatrix`、张量工具 |
| `defaults.jl` | 全部默认超参数（`Defaults`） |
| `mpohamiltonian/` | `SchurMPOTensor`、`MPOHamiltonian`、长程衰减项、`timeevompo`（WI/WII/ComplexStepper） |
| `adt/` | ADT（MPS）类型、正交化、SVD 压缩、`mult!`、积分/交换门 |
| `pt/` | `ProcessTensor`（MPO）类型、正交化、乘法 |
| `adtlattices/` | 实/虚/混合时间 ADT 格点与索引映射、Fock 排序 |
| `ptlattices/` | 实/虚/混合时间 PT 格点、`integrate`（PT 求值） |
| `correlationfunction.jl` | 浴关联函数包装（`IndexCorrelationFunction`、分支关联） |
| `influencefunctional/` | 影响泛函构造：`PartialIF`（bond-dim-2）、`TranslationInvariantIF`（XTRG 风格）、`TDVPIF`（二阶单格点 TDVP 虚时间流）、`PTPartialIF`、PT 的 TTI-IF |
| `tdinfluencefunctional/` | 含时耦合（`*TdHyb`）影响泛函 |
| `boundarycondition.jl` | 初态/边界条件 `boundarycondition!`、`initialstate!` |
| `models/` | 系统哈密顿量/Lindblad 算符、`sysdynamics` |
| `observables/` | 环境缓存、期望值、配分函数、转移矩阵 |

## 2. 核心数据结构

### 2.1 ADT（Augmented Density Tensor / MPS）

定义在 `src/adt/def.jl`：

```julia
struct ADT{T<:Number, R<:Real} <: Dense1DTN{T}
    data::Vector{Array{T, 3}}            # 位点三阶张量
    s::Vector{Union{Missing, Vector{R}}} # Schmidt 向量（奇异向量）
    scaling::Ref{Float64}                # 全局标度因子
end
```

- 张量约定（三腿）：腿 1 = 左辅助（bond）指标，最左位点维数为 1；腿 2 = 物理指标（1 对应 $|0\rangle$，2 对应 $|1\rangle$）；腿 3 = 右辅助指标，最右位点维数为 1。
- `s` 为每个 bond 的奇异值向量（长度 = 位点数 + 1），用于 `iscanonical` 校验与 `easy_swap!` 交换门。
- `scaling` 记录整体标度（见 §2.3）。

### 2.2 ProcessTensor（MPO）

`ProcessTensor` 为四腿张量列（`Vector{Array{T,4}}`），腿顺序为 `(左bond, 物理in, 物理out, 右bond)`，表示影响泛函算符链（MPO）。`ADT` 可看作 MPO 的特例（out 腿缩并到 in 腿，即"对角化"状态）。

### 2.3 标度因子（scaling）机制

为了数值稳定性，所有正交化/截断步骤都避免直接归一化张量，而是把范数累积到全局 `scaling` 中（`src/adt/orth.jl`）：

- `_renormalize!(psi, r, normalize)`：计算 `r` 的范数 `nr`，若 `!normalize` 则调用 `_rescaling!` 把 `nr^(1/L)` 乘入 `scaling`，再 `lmul!(1/nr, r)`；
- `_rescaling!(psi, n)`：`scaling *= n^(1/L)`。

这样 bond 张量保持 O(1) 量级，避免指数增长/衰减导致的上溢/下溢。`mult!` 的结果标度为 `scaling(x) * scaling(y)`。

## 3. 轮廓与格点

### 3.1 时间轮廓与分支

- 虚时间（Matsubara）：分支 `:τ`，区间 $[0, \beta]$，标量类型 `Float64`；
- 实时间（Keldysh）：分支 `:+` 与 `:-`，区间 $[0, t]$，标量类型 `ComplexF64`；
- 混合（Kadanoff-Baym）：虚时间分支 + 实时间上下两支。

`branch`/`index` 由 `QuAPI` 提供并导入（`src/TEMPO.jl`）。每个格点还关联 `TimeOrderingStyle`（`TimeAscending`/`TimeDscending`）与 `LayoutStyle`（当前只有 `TimeLocalLayout`，即每个时间步的状态局部排列，便于施加时间局域算符）。

### 3.2 Fock 排序（FockOrdering）

Fock 排序定义了格点上 Grassmann 算符的排列约定（`src/adtlattices/fockordering.jl`）：

| 类型 | 轮廓 | 含义 |
|---|---|---|
| `M2M1`（别名 `MM`） | 虚时间 | 每个虚时间步含两个序 M2、M1；虚支时间降序 |
| `M2m2M1m1`（别名 `MmMm`） | 实时间 | 每个实时间步含两个序及其共轭 m2、m1；实支时间降序 |
| `M2M1_m1M1m2M2` | 混合 | 实支时间升序，虚支时间降序 |

### 3.3 格点类型与索引映射

实时间 ADT（`src/adtlattices/realtime.jl`）：

```julia
struct RealADTLattice1Order{O<:RealFockOrdering} <: RealADTLattice{O}
    δt::Float64; d::Int; N::Int; ordering::O
end
```

派生属性：`t = N*δt`、`Nt = N`、`k = kt = N+1`、`length = 2k`。索引映射（降序，`+1` 为边界 Grassmann 数留位）：

```julia
index(i, branch=:+) = TL - 2i + 1    # TL = 2(N+1)
index(i, branch=:-) = TL - 2i + 2
```

即物理位点按 `(N,+) (N,−) (N−1,+) (N−1,−) … (1,+) (1,−)` 排列，位置 1、2 留给初态边界。

虚时间 ADT：`length = k = N+1`，`index(i) = k+1-i`（降序），最右端（i=1）为真空边界。

PT 格点（`src/ptlattices/`）：
- `ImagPTLattice1Order`：`length = N`，`index(i) = N-i+1`，仅 `:τ` 分支，标量类型 `Float64`；
- `RealPTLattice1Order`：`length = 2N`，`index(i, :+) = 2(N-i)+1`、`index(i, :-) = 2(N-i)+2`（降序，`index(N, :+)=1`，`index(1, :-)=2N`）；
- `MixedPTLattice1Order`：`Nt` 个实时间步 + `Nτ` 个虚时间步，`index` 在实/虚支之间交错。

`indexmappings(lattice)` 返回 Dict{(timestep, branch) => global_index}，`vacuumstate(T, lattice)` 构造全 1 的真空 ADT。

## 4. 浴关联函数与指数展开

`correlationfunction.jl` 把关联函数包装为 `IndexCorrelationFunction`（内含 `CorrelationMatrix` 的 ηᵢⱼ 矩阵，由 `exponential_expansion` 生成），实时间用 `branch(corr, :+, :-)` 等取四个符号化分量 η⁺⁺、η⁺⁻、η⁻⁺、η⁻⁻。

指数展开（`ExpExp` 包，经 TEMPO 重导出）：
- `OverDeterminedProny`：将关联函数展开为有限项指数和（最小二乘 Prony 方法）；
- `DeterminedProny`：确定型 Prony 方法变体；
- `expand_decayterm(decayterm, alg=...)` 返回 `(η₁, η₂, …)` 系数列；
- `expansion_error` 估计展开误差。

## 5. 对角耦合影响泛函（PartialIF / ADT 路径）

### 5.1 `partialif_densemps`：bond 维数 2 的部分影响泛函

算法来自 Strathearn et al. (2018)，实现在 `src/influencefunctional/partialif/util.jl`：

- 输入：行位置 `row`、列位置 `cols`、对角算符 `op`（`z`）、系数 `coefs`（η 值）；
- 构造单行 MPS（bond 维数 = `d`，对角耦合 d=2 时为 2）：
  - 行位点处放置"单体"张量 `exp(ηₛₛ z²)`（对角矩阵，`tmp[i,i,i]`）；
  - 其余位点放置"两体"门 `exp(ηᵢⱼ z⊗z)`（`tmp[i,:,i] = m[i,:]`，bond 指标承载物理态）；
  - 首/末位点为边界张量 `(1,d,d)`/`(d,d,1)`；
- `_fit_to_full` 把 MPS 嵌入整条格点：行前后补全 1 的真空张量 `ones(1,d,1)`，行区间内未命中的位置放 `Σᵢ δ` 型恒等张量（`tmp[:,i,:] = I`）。

### 5.2 `partialif_naive`

朴素版本：对每个列索引依次 `apply!` 指数门 `exp(coef * op⊗op)` 后 `canonicalize!` 截断，最终得到部分 IF（同文件）。`coef` 取实时间用 ηᵢⱼ（含符号），虚时间用 ηᵢⱼ/2 对半分配。

### 5.3 `hybriddynamics!`

（`src/influencefunctional/partialif/realtime.jl` 及虚/混合时间对应实现）逐行构造部分 IF 并依次乘入：

```julia
for i in 1:Nt, b1 in branches(lattice)
    pos1 = index(lattice, i, branch=b1)
    # 收集该行与所有列的系数
    tmp = partialif_densemps(ds, pos1, pos2s, op, coefs)
    mult!(gmps, tmp, trunc=trunc)
end
```

每个 `(i, b1)` 行的部分 IF 与全局 `gmps` 相乘并 SVD 压缩（`DefaultITruncation`）。这是文献 Eq. (2)-(8) 的 ADT 实现。

## 6. 平移不变影响泛函（TTI-IF，XTRG 风格）

对应文献附录中"平移不变 + 指数展开 + MPO 时间演化"方案，入口为 `influenceoperator`/`influenceoperatorexponential`/`differentialinfluencefunctional`（`src/influencefunctional/ttiif/`）。

### 6.1 指数展开与 `SchurMPOTensor`

`adt_ti_mpotensor`（`src/influencefunctional/ttiif/adt/imag.jl`）：

```julia
m1 = GenericDecayTerm(op1, op2, corr.ηⱼₖ[2:end])   # 长程（跨步）耦合
m2 = GenericDecayTerm(op2, op1, corr.ηₖⱼ[2:end])
m1s = expand_decayterm(m1, alg)                     # Prony 展开
m2s = expand_decayterm(m2, alg)
h1  = (corr.ηₖⱼ[1] + corr.ηⱼₖ[1]) .* (op1 * op2)    # 同时间（对角块）项
return SchurMPOTensor(h1, vcat(m1s, m2s))
```

`SchurMPOTensor` 为块三角算符结构（`src/mpohamiltonian/schurmpo/`）：`D = h1`（对角块，同时间耦合），`A = {decay terms}`（对角块，长程耦合），`B`、`C` 为上下三角连接块。

### 6.2 时间演化：WI / WII / ComplexStepper

`timeevompo(m, dt, alg)`（`src/mpohamiltonian/schurmpo/w1w2.jl`），方案来自 Zaletel et al., arXiv:1407.1832：

- **WI**（一阶）：`WD = I + dt·D`，`WB = B·√δt`，`WC = C·√δt`（`_sqrt2` 对负 dt 给出 `(√|dt|, −√|dt|)`），再组装为稀疏 MPO 张量；
- **WII**（一阶，精度更高）：构造 `4d×4d` 块矩阵

  ```
  [[Ddt,  0,   0,   0 ],
   [√δ₂C, Ddt, 0,   0 ],
   [√δ₁B, 0,   Ddt, 0 ],
   [A,    √δ₁B,√δ₂C,Ddt]]
  ```

  取 `exp` 后按列提取块 `WC`、`WB`、`WA`，`WD = exp(Ddt)`；
- **ComplexStepper**：`dt₁ = (1−i)dt/2`、`dt₂ = (1+i)dt/2`，两次一阶步复合为二阶精度，返回演化前/后的两组 MPO（`timeevompo` 返回 `(U₁, U₂)`）。

### 6.3 映射到格点：`_fit_to_lattice`

虚时间 ADT 中，`_fit_to_lattice`（`src/influencefunctional/ttiif/adt/imag.jl`）把 3 个 MPO 张量按平移不变模式铺到 `L = N+1` 个位点：

- 位置 1（j=N+1）← `mpstensors[1]`（左边界张量）；
- 位置 2…N−1（j=N−1…3）← `mpstensors[2]`（体张量，重复）；
- 位置 N（j=2）← `mpstensors[3]`（右边界张量）；
- 位置 N+1（j=1）← `ones(1, d, 1)`（真空）。

`_tompsj` 把 MPO 张量的 out 腿用全 1 向量缩并（`mpoj[1,2,3,4]*a[4]`），将算符化为"态"张量以便嵌入 ADT。实时间 ADT 类似但返回 4 个分支 MPO 的元组。

### 6.4 实时间分支 MPO 与微分影响泛函

PT 实时间 `influenceoperator`（`src/influencefunctional/ttiif/pt/real.jl`）返回 4 个分支 MPO `(η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻)`；`influenceoperatorexponential` 对每个分支先 `timeevompo` 再拟合，`FirstOrderStepper` 返回 4 个、`ComplexStepper` 返回 8 个。`differentialinfluencefunctional` 依次相乘（PT 分支顺序 `h2*h1`、`h3*…`、`h4*…`）得到单个时间步的完整微分影响泛函。

## 7. 非对角耦合 PT 框架

对应文献"有效哈密顿量 + 拷贝张量"的非对角耦合方案。

### 7.1 有效哈密顿量 H_eff

系统 + 浴的非对角耦合（如 JC 模型）通过引入辅助"拷贝"自由度化为对角形式：`H_eff` 中拷贝系统态与浴算符配对，使影响泛函仍为指数门乘积。工具包用 PT（MPO）表示这条路径的影响泛函。

### 7.2 `pt_ti_mpotensor` 与符号化关联函数

对每个分支组合构造 MPO 张量：

```julia
mpoj1 = pt_ti_mpotensor(η⁺⁺, op1, op2, :+, :+, algexpan)
mpoj2 = pt_ti_mpotensor(η⁺⁻, fused_op(op1, :+), fused_op(op2, :-), :+, :-, algexpan)
mpoj3 = pt_ti_mpotensor(η⁻⁺, fused_op(op1, :-), fused_op(op2, :+), :-, :+, algexpan)
mpoj4 = pt_ti_mpotensor(η⁻⁻, transpose(op1), transpose(op2), :-, :-, algexpan)
```

- 对角分支（++, −−）直接用原算符（−− 分支取转置，对应共轭传播子）；
- 非对角分支（+−, −+）用 `fused_op` 构造"融合"两腿算符：`:+` 支作用 `op1`、identity 作用在 − 支；`:-` 支作用 identity、`transpose(op1)` 作用在 − 支。

### 7.3 `fused_op` 与 `split_mpotensor`

- `fused_op(op, f)`（`src/influencefunctional/ttiif/pt/real.jl`）：用恒等嵌入张量 `f = reshape(I_{d²}, d², d, d)` 把单腿算符提升为双腿 `(d²×d²)` 算符；
- `split_mpotensor(mpoj, trunc)`：把 4 腿 MPO 张量按物理维 `d²` 重构为 6 腿张量，SVD 拆成 `u`、`v` 两块（`u` 带 `√s`，`v` 带 `√s`），用于非对角分支在两条支路的不同位置放置耦合。

### 7.4 对角/非对角分支的格点拟合

- `_fit_to_lattice_diag`：对角分支只在一个支路位置放 MPO 张量，另一支路位置放 `vd2 ⊗ I₂` 恒等张量；
- `_fit_to_lattice_offdiag`：非对角分支在 `pos1`（=min 位置）放 `u` 块、`pos2` 放 `v` 块，中间位置填恒等张量（`band_boundary` 给出每个时间步 ± 分支的位置区间）。

由此，非对角耦合的指数门被"分裂"到上、下两支的格点位置上，物理上等价于有效哈密顿量路径的拷贝张量。

### 7.5 PT→ADT 转换

`conversions.jl` 中的 `copytensor` 构造拷贝张量 `m[i,i,i]=1`；PT→ADT 的完整转换（`toadt`）在当前版本中处于注释/占位状态（文献方法已有实现线索，见文件头部注释）。

## 8. 系统动力学与边界条件

- `boundarycondition!`（`src/boundarycondition.jl`）：在格点边界上施加初态（`initialstate!`）。实时间轮廓在位置 1、2 放置初始密度矩阵的 Grassmann 表示；虚时间在两端放置真空/热态边界；
- `sysdynamics`/`sysdynamics!`（`src/models/`）：把系统的时间演化算符（一阶指数门 `exp(-iH δt)`，或 Lindblad 传播子）乘入 ADT/PT，对应文献 Eq. (9)-(10) 的系统传播子；
- `ImpurityHamiltonian`/`ImpurityLindbladian`：由 `ImpurityModelBase` 定义系统哈密顿量/耗散算符，`sysdynamics` 据其生成传播子并就地作用。

## 9. MPO 哈密顿量机制

`MPOHamiltonian`（`src/mpohamiltonian/mpohamiltonian.jl`）为一列 `SchurMPOTensor`，用于：

- 构造系统哈密顿量的 MPO 表示（`tompotensors`）；
- 时间演化 `timeevompo(m, dt; alg=WII())`（逐位点演化）；
- 长程耦合项：`GenericDecayTerm`（通用指数衰减）、`PowerlawDecayTerm`（幂律，如 XXZ 长程）。

`get_A/B/C/D` 提取 `SchurMPOTensor` 的块结构，`_SiteW_impl` 组装为 `SparseMPOTensor`。

## 10. 张量网络算法

### 10.1 正交化

（`src/adt/orth.jl`）`Orthogonalize` 配置 `orth`（QR/SVD）、`trunc`、`normalize`、`verbosity`：

- `leftorth!`：从左到右 `tqr!`/`tsvd!`，R/S 归一化后吸收进下一位点；
- `rightorth!`：从右到左 `tlq!`/`tsvd!`，L 归一化后吸收进前一位点；
- `canonicalize!`：先用 QR 左正交（不截断），再用 `alg` 指定的 SVD 右正交并截断。注释明确警告：对 ADT/ProcessTensor **不要启用 normalize**（重归一化会破坏影响泛函的标度语义）。

QR 路径中 truncation 无效（`@warn`）。截断误差按相对误差 `sqrt(ε²/(n²+ε²))` 报告。

### 10.2 MPS 乘法 `mult!` / `mult`

（`src/adt/mult/svdmult.jl`）`mult!(x, y)` 算法：

1. 左端 `tqr!`（`tie(n_fuse(...), ...)` 合并指标）；
2. 从左到右迭代：`tmp = r ⊗ x[i] ⊗ y[i]` 收缩 → `n_fuse` → `tqr!` → `_renormalize!`；
3. 右端收尾后 `_rightorth!(x, SVD(), trunc)` 从右向左 SVD 截断；
4. `setscaling!(x, scaling(x)*scaling(y))`。

`mult(x, y) = mult!(copy(x), y)` 为不修改输入的版本。`DMRGMult1`（`DMRGAlgorithm`）提供带 `initguess`（默认 `:svd`）的变体，配合 `SVDCompression` 的 `D`/`tol` 截断。

### 10.3 截断方案

- `truncdim(D)`：仅按最大 bond 维截断；
- `trunccutoff(ε)`：按累积截断误差截断；
- `truncdimcutoff(; D, ε)`：两者结合；`NoTruncation()` 不截断。

### 10.4 转移矩阵

`ADTTransferMatrix`（`src/observables/adt/transfer.jl`）把一组 ADT 组装为转移算符：

- `left * m`：从左到右逐位点 `transfer_left` 收缩（`lmul!(scaling, …)` 计入标度）；
- `m * right`：从右到左；
- `TransferMatrix(states...)` / `TransferMatrix(j, states...)`：整体或单点版本（单点算符用）。

## 11. 观测量计算

### 11.1 环境缓存

`environments(lattice, A, B...)`（`src/observables/`）构造环境缓存：

- `ADTExpectationCache{A, Bs, lattice, hleft, hright}`：`hleft[i]` 为前 i 个位点的左环境（从左到右 `left * TransferMatrix(i, As...)` 累积），`hright[i]` 为右环境；`Zvalue = only(hleft[end])`；
- `PTExpectationCache`：PT 版本，实时间轮廓接受关键字 `ρ₀`（初态密度矩阵）作为右边界条件。

### 11.2 期望值

- `expectation(m, cache)`：对 `copy(cache.A)` 施加算符 `m`（`apply!`），在算符支撑区间上从右向左收缩转移矩阵，最后 `contract_center(left, right)` 缩并两端环境；
- `expectationvalue(m, cache) = expectation(m, cache) / Zvalue(cache)`：归一化期望值 $\langle m\rangle = \mathrm{tr}_{\text{bath}}[\bar\psi m \psi]/Z$；
- `Zvalue(cache)`：配分函数 Z。

## 12. 误差来源与超参数

文献中五类误差在代码中的对应：

| 误差来源 | 超参数 | 说明 |
|---|---|---|
| 时间离散 | `δt`（`δτ`） | 一阶 Trotter/影响泛函离散误差 |
| SVD 截断 | `truncdim` / `trunccutoff` / `truncdimcutoff` | bond 维数 `χ` |
| 指数展开 | `algexpan` | Prony 展开项数与容差（默认 `OverDeterminedProny(n=15, tol=1e-4)`） |
| 平移不变细化 | `k`（默认 5）、`fast` | `fast=true`：先构造宽度 `dt/2^k` 的微分 IF，再树二分平方 k 次得到全长影响泛函 |
| 系统传播子 | `algevo`（`WII()`）、`algmult`（`DefaultMultAlg`） | MPO 时间演化与乘法压缩精度 |

默认值见 `src/defaults.jl`：`DefaultTruncation = truncdimcutoff(D=100, ϵ=1e-14)`、`DefaultITruncation = truncdimcutoff(D=200, ϵ=1e-10)`、`DefaultKTruncation = truncdimcutoff(D=1000, ϵ=1e-10)`、`DefaultIntegrationTruncation = DefaultMPOTruncation = truncdimcutoff(D=10000, ϵ=1e-12)`；`TranslationInvariantIF(; algexpan=OverDeterminedProny(n=15, tol=1e-4), algevo=WII(), algmult=DefaultMultAlg, k=5, fast=true)`。

## 13. 数据流总览

计算一条完整的时间演化动力学：

```
corr = bath 关联函数（ImpurityModelBase 谱函数/相关函数）
  └─ exponential_expansion → ηᵢⱼ 系数
lattice = ADTLattice/PTLattice（轮廓 + Fock 排序 + 索引映射）
hyb = AdditiveHyb（对角）或 NonDiagonalHyb/GeneralHybStyle（非对角）
  ├─ 对角：PartialIF（bond dim 2）逐行乘法，或 TranslationInvariantIF（MPO 指数展开 + WI/WII 演化）
  └─ 非对角：PT 路径，4 分支 MPO（η⁺⁺…η⁻⁻）+ fused_op/split_mpotensor + 分支相乘
sysdynamics!：系统传播子（H 或 Lindblad）乘入
boundarycondition!：初态
observables：environments 缓存 → expectationvalue/expectation/Zvalue
```

核心循环的标度关系：ADT/PT 的位点数 `L` 与时间步数线性增长，bond 维数受截断方案控制；TTI 方案把单步微分影响泛函（MPO）按平移不变模式复制到所有位点，大幅降低影响泛函构造成本。
