# 接口调整说明（2026-09-17）：截断方案与 MPS 算法接口统一（SVDCompression/DMRG1 参数化、默认截断收敛）

对齐主流库（MPSKit / ITensor / TeNPy）的做法，把压缩算法的截断参数统一为 `TruncationScheme` 对象，并收敛默认截断常量。全套测试通过（0 Fail / 0 Error）。

## 截断方案（`src/tensorops/truncation.jl`）

| 旧名（已删除） | 新名 | 说明 |
|---|---|---|
| `TruncationDimCutoff` | `TruncateDimCutoff` | 类型更名，与其余 `Truncate*` 命名对齐；便利构造函数 `truncdimcutoff(D, ϵ[, add_back])` **不变** |

- `trunccutoff` 新增位置参数构造 `trunccutoff(ϵ::Real)`（与关键字形式 `trunccutoff(; ϵ)` 等价）。

## MPS 算法（`src/algorithms.jl`）

### `SVDCompression`：`D`/`tol` 字段 → 参数化 `trunc` 字段

| 旧接口 | 新接口 |
|---|---|
| `SVDCompression(; D=Defaults.D, tol=Defaults.tol, verbosity=0)` | `SVDCompression(; trunc=truncdimcutoff(D=Defaults.D, ϵ=Defaults.tol, add_back=0), verbosity=0)` |
| `SVDCompression(trunc::TruncationDimCutoff; verbosity=0)` | `SVDCompression(trunc::TruncationScheme; verbosity=0)`（接受**任意** `TruncationScheme`） |
| `alg.trunc`（getproperty 合成）/ `get_trunc(alg)` / `alg.ϵ` | `alg.trunc`（真实字段） |

- 结构体变为 `SVDCompression{T<:TruncationScheme}`；`Base.similar(; trunc, verbosity)` 同步参数化。

### `DMRG1`：`trunc` 参数化但必须携带键维 `D`

| 旧接口 | 新接口 |
|---|---|
| `DMRG1(trunc::TruncationDimCutoff; ...)` | `DMRG1(trunc::TruncationWithD; ...)`，其中 **`TruncationWithD = Union{TruncateDim, TruncateDimCutoff}`** |
| `DMRG1(; trunc::TruncationDimCutoff=DefaultITruncation, ...)` | `DMRG1(; trunc::TruncationWithD=DefaultITruncation, ...)` |

- 原因：`iterativemult` 的初始猜测（`:svd`/`:rand`/`:pre`）需要 `D` 信息，实现改用 `alg.trunc.D`。
- **删除** `Base.getproperty(::DMRGAlgorithm, :D/:ϵ)` 访问器（对无 `D`/无 `ϵ` 的方案无定义）；`Base.similar` 同步参数化。

## 默认截断常量收敛（`src/defaults.jl`）

| 常量 | 变更 |
|---|---|
| `DefaultKTruncation` | `truncdimcutoff(D=1000, ϵ=1e-10)` → **`trunccutoff(Defaults.tolgauge)`** |
| `DefaultIntegrationTruncation` | **已删除**，原用点改用 `DefaultKTruncation` |
| `DefaultMPOTruncation` | **已删除**，原用点改用 `DefaultKTruncation` |
| `DefaultTruncation` / `DefaultITruncation` / `DefaultMultAlg` | 不变 |

替换位置：`boundarycondition!`、`_permute!`（ADT/PT linalg）、TTIIF 的 `_fit_to_lattice_diag/_offdiag`（adt/pt real）、TDVPIF 的 H 压缩（`_tdvpif_hamiltonian`）。

## 其他删除

- `src/observables/correlations.jl`：整个文件删除（从未被 include 的死代码，`correlation` 函数无处分发使用；两点关联测量请用 `ADTTerm` 多点形式 + `apply!`/`integrate`，见 manual §Observables）。`docs/src/api.md`、`docs/src/manual.md` 同步清理。

## 迁移指南

- `SVDCompression(D=χ)` → `SVDCompression(truncdimcutoff(D=χ, ϵ=Defaults.tol))`（或按需写 `truncdimcutoff(D=χ, ϵ=...)`）。
- `SVDCompression(D=χ, tol=ε)` → `SVDCompression(truncdimcutoff(D=χ, ϵ=ε))`。
- `DMRG1(trunc)` / `DMRG1(trunc=...)`：`trunc` 需为 `truncdim(D)` 或 `truncdimcutoff(D, ϵ)`。
- `alg.D` / `alg.ϵ` 访问改为 `alg.trunc.D` / `alg.trunc.ϵ`。
- `DefaultIntegrationTruncation` / `DefaultMPOTruncation` → `DefaultKTruncation`。
- `TruncationDimCutoff` → `TruncateDimCutoff`。

## 测试

- `test/api/truncation.jl` 新增 `SVDCompression / DMRG1` testset：参数化构造（全部 `TruncationScheme` / 仅含 `D` 方案）、`MethodError` 拒绝、`similar` 保持/覆盖、默认值；`test/api/mps.jl` 的 multiplications 增加 `truncdim`/`trunccutoff` 方案下的乘法正确性；`test/api/influenceoperator.jl` 的 `SVDCompression(D=50)` 调用点迁移。
- 全套测试通过（0 Fail / 0 Error）。

---

# 接口调整说明（2026-09-16）：迭代乘法（DMRG1）收敛判据与 sweep!/iterative_compute! 接口统一

参考 MPSKit / ITensor / TeNPy / quimb / block2 的主流做法，统一迭代乘法（`mult(x, y, alg::DMRGAlgorithm)`，ADT 与 PT 共用）的 loss 度量与收敛判据接口。全套测试通过（1093/1093，0 Fail / 0 Error）。

## loss 的定义与实现细节

- 迭代乘法求解变分问题 `z ≈ w ≡ x·y`（未截断乘积网络），loss 为 `F(z) = ‖w − z‖²`；single-site ALS 逐站点最优更新保证 F 单调不增（数值验证至机器精度）。
- 每个 sweep 中，站点 `j` 的局部最优张量 `mpsj_j = L·w_j·R`（Riesz 代表元），其范数 `‖mpsj_j‖` 即该站点的 loss 值；定点处满足 `‖mpsj_j‖² = ‖w‖² − F*`，故残差衡量“已捕获质量”：sweep 内严格单调上升至常数 `√(‖w‖²−F*)`，loss 相应单调下降。
- 注意：`sweep` 过程中 QR/LQ 的尾巴（中心矩阵）被丢弃，存储张量存在等距类内的规范漂移；因此基于存储张量差的判据（如 MPSKit `approximate!` 的 `‖AC′−AC‖/‖AC′‖`）不适用于本实现，已验证其恒 ≈ 1 不收敛。

## 接口变更

| 函数 | 旧接口 | 新接口 |
|---|---|---|
| `sweep!(m, alg)` | 返回该 sweep 所有逐站点 loss 值 `‖mpsj_j‖` 的 vector（行为不变，此处明确为接口约定） | 不变 |
| `iterative_compute!(m, alg)` | 返回逐 sweep 的判据值向量 | **返回逐 sweep 的末位 loss 值向量** `kvals[t] = 第 t 轮 sweep 的最后一个残差`（该序列单调上升至定点 `√(‖w‖²−F*)`） |
| 收敛判据 | 旧：sweep 内残差的相对涨落 `std/mean`（本轮之前）→ 相邻 sweep 残差逐站点相对差 max（本轮中间版） | **相邻 sweep 的末位 loss 相对差** `δ = |r_last(t) − r_last(t−1)| / max(r_last(t), r_last(t−1)) < tol`；首轮强制运行（`delta = 2*tol`） |
| `iterative_error_2` | 旧判据工具（sweep 内残差 `std/mean`） | **已删除**（无调用点） |

## 文档

- `docs/src/internals.md` 新增 §10.3（迭代乘法）：loss 泛函、残差含义与单调性、sweep 结构（QRpos/LQpos、finalize 截断、初始 guess）、收敛判据及与 MPSKit 判据不兼容的原因；原 §10.3/§10.4 顺延为 §10.4/§10.5。

## 测试

- 全套测试通过（1093/1093，0 Fail / 0 Error）；debug 脚本 `debug/mult/run_mult_debug.jl`（不入库）同步新判据并复验 72 算例：收敛轮数与原判据同量级、物理结果逐数一致。

---

# 接口调整说明（2026-09-09）：TTIIF 影响算子函数更名与虚时返回值统一为元组

对齐 GTEMPO 的命名。全套测试通过（0 Fail / 0 Error）。

## 函数更名

| 旧名（已删除） | 新名 | 说明 |
|---|---|---|
| `influenceoperator` | `influenceoperators` | 分支影响算子 MPO 组 |
| `influenceoperatorexponential` | `influenceoperatorsteppers` | 单步演化后的影响算子（stepper）组 |
| `differentialinfluencefunctional` | `influenceoperatorstepper` | 单步差分影响泛函（各分支 stepper 的乘积） |

## 返回值统一

- **虚时 `influenceoperators` 的返回值由裸 MPO 改为 1 元组 `(mpo,)`**，与实时的 4 元组约定一致；调用方使用 `only(...)` 或 `mpo, = ...` 解包（TDVPIF 的 `hybriddynamics!` 内部已同步改为 `only(influenceoperators(...))`）。
- `influenceoperatorsteppers` 的返回值本就是元组（虚时 FirstOrder 1 元组 / ComplexStepper 2 元组；实时 4/8 元组），不变。

## 迁移指南

- `influenceoperator(...)` → `only(influenceoperators(...))` 或解包；实时多返回值调用无需改动解包方式，仅函数名变化。
- `influenceoperatorexponential` / `differentialinfluencefunctional` 直接改名为 `influenceoperatorsteppers` / `influenceoperatorstepper`，参数不变。

---

# 接口调整说明（2026-09-09）：含时杂质哈密顿量支持

参考 GTEMPO 的 `QuenchedImpurityHamiltonian` / `TdImpurityHamiltonian`，为 TEMPO 增加含时杂质哈密顿量支持（矩阵约定）。全套测试通过（含新增 testset，0 Fail / 0 Error）。

## 新增类型（`src/models/def.jl`，均已导出）

- **`AbstractImpurityHamiltonian <: AbstractImpurityOperator`**：单元型杂质哈密顿量模型的抽象父类型；`ImpurityHamiltonian` 改为继承它（对外接口不变）。
- **`QuenchedImpurityHamiltonian(hτ, ht)`**（quench 协议）：虚时间分支（`:τ`）以 `hτ` 演化，实时间分支（`:+`/`:-`）以 `ht` 演化；两矩阵尺寸必须一致。
- **`TdImpurityOp(m, f)`**：含时项，`t` 时刻贡献 `f(t)·m`。
- **`TdImpurityHamiltonian(hτ, ht, [TdImpurityOp...] )`**：虚时间分支以常数 `hτ` 演化；实时间分支以逐步哈密顿量 `ht + Σₖ fₖ(t)·mₖ` 演化。调用 `model(t)` 返回 `t` 时刻实分支的哈密顿量矩阵。

## 逐步传播子接口

- 新增 `propagator(model, lattice, branch, j, N)`：第 `j`/`N` 步的传播子；默认回退到（与步无关的）`propagator(model, lattice, branch)`，`TdImpurityHamiltonian` 覆盖之。步时间约定与 GTEMPO 一致：forward 分支 `t = (j-1)δt`，backward 分支 `t = (N-j)δt`。
- `TdImpurityHamiltonian` 的实分支传播子逐步重算（`exp(-im·δt·H(t))`）；虚分支/常数模型仍为分支常数传播子。

## dynamics 改动

- ADT / PT 的 `sysdynamics!`、`sysdynamics_forward!/backward!/imaginary!` 的模型签名放宽为 `AbstractImpurityHamiltonian`；`ImpurityHamiltonian` 仍走原有的"传播子提升到循环外"的常数值路径（行为不变）。
- 新增通用逐步循环 `_sysdynamics_util!`（ADT 逐 gate `apply!` + 每步 canonicalize；PT 逐 `ContourOperator` `apply!` + 末尾 canonicalize），支持 `Quenched`（各分支内为常数）与 `Td` 模型。
- 含时模型三种轮廓（imag / real / mixed）均支持；暂不支持与 `ContourOperator` 插入（`sysdynamics!` 的 cts 变体）及 Lindblad 耗散组合。

## 测试

- 新增 `test/models/tdimpurity.jl`（并入 `test/models/models.jl`）：quench 各分支等价于对应常数模型（ADT/PT × imag/real/±branch）；`Td` 虚分支 == `hτ` 模型；常数 `TdImpurityOp` == `ht + c·m` 常数模型（ADT/PT real + mixed 轮廓路由）；真实含时演化与测试内独立构造的逐步精确传播子逐一比对（固定步时间约定）；`model(t)` 数值校验。

---

# 接口调整说明（2026-09-09）：DMRG 迭代乘法算法更名与算法定义归位

本轮对齐 GTEMPO 的算法类型层级：删除 `DMRGMultAlgorithm`，`DMRGMult1` 更名为 `DMRG1`，并把算法类型定义集中到 `src/algorithms.jl`、默认值集中到 `src/defaults.jl`。行为不变，全套测试通过（93 个 testset，0 Fail / 0 Error）。

## 类型层级

| 旧名（已删除） | 新名 | 说明 |
|---|---|---|
| `DMRGMultAlgorithm` | —（删除） | 原 `DMRGMult1` 的父类型；相关方法（`mult`/`mult!`/`iterativemult`/`compute!`/`sweep!`/`finalize!`）改分派到 `DMRGAlgorithm` |
| `DMRGMult1` | `DMRG1` | 单点 DMRG 迭代乘法；直接继承 `DMRGAlgorithm`。字段与构造器不变（`trunc`, `maxiter`, `tol`, `initguess ∈ {:svd, :pre, :rand}`, `verbosity`, `callback`） |

注：TEMPO 中不存在 `DMRGMult2`，故本轮仅涉及 `DMRG1`。

## 定义位置迁移

- **迁入 `src/algorithms.jl`**（与 `MPSAlgorithm`/`DMRGAlgorithm`/`SVDCompression` 同处）：
  - `MatrixProductOrthogonalAlgorithm` 抽象类型与 `Orthogonalize{A<:Union{QR, SVD}, T<:TruncationScheme}` 及其三个构造器（原 `src/adt/orth.jl`）；
  - `AllowedInitGuesses` 常量、`DMRG1` 结构体/构造器/`Base.similar`（原 `src/adt/mult/iterativemult.jl`）；
  - `Base.getproperty(::DMRGAlgorithm, :D/:ϵ)`（原分派于 `DMRGMultAlgorithm`）。
- **迁入 `src/defaults.jl`**：`DefaultMultAlg = DMRG1(DefaultITruncation)`（原 `src/adt/mult/mult.jl`；置于 `DefaultITruncation` 之后，无初始化顺序问题）。

## 迁移指南

- `DMRGMult1(...)` → `DMRG1(...)`；`DMRG1` 已导出，`DMRGMult1` 不再导出。
- 类型注解 `::DMRGMultAlgorithm` → `::DMRGAlgorithm`。
- `Orthogonalize`、`leftorth!`/`rightorth!`/`canonicalize!`、`SVDCompression` 等接口不变。

---

# 接口调整说明（2026-09-03 / 09-04）

这一波改动重构了底层张量分解模块的文件组织与函数接口，并统一了原地 / 非原地版本的语义。

## 文件结构

| 原路径 | 新路径 | 说明 |
|---|---|---|
| `src/auxiliary/` | `src/tensorops/` | 目录重命名 |
| `src/auxiliary/tensorops.jl` | `src/tensorops/tensorfactorizations.jl` | 张量分解函数（`tsvd`/`leftorth`/`rightorth` 等） |
| `src/auxiliary/auxiliary.jl` | `src/tensorops/tensorops.jl` | 聚合 include 入口 |
| `src/auxiliary/mpsalgs.jl` | `src/algorithms.jl` | MPS/DMRG 算法（`MPSAlgorithm`、`SVDCompression`），移出至 `src/` 顶层 |

同步更新：`src/TEMPO.jl` 的 include 与 export、`docs/src/api.md` 的 `@autodocs` Pages 路径、`manual.md`/`internals.md` 的目录说明。

## 函数接口变更

### `tsvd!` / `tsvd`

- **删除 `stable_tsvd!`**。
- `tsvd!(a; trunc, alg)`（矩阵版）改用 MatrixAlgebraKit 的 `svd_compact!`，是真正的**破坏性**版本：输入矩阵用作 workspace，可能被覆写（原实现基于 out-of-place 的 `svd_compact`，名不副实）。
- 新增 `alg` 关键字选择 SVD 驱动：
  - `SDD()`（默认）：divide-and-conquer（`SafeDivideAndConquer` / LAPACK `gesdvd`，自带 QR-iteration 回退，即原 `stable_tsvd!` 的稳健行为）；
  - `SVD()`：QR iteration（`QRIteration` / LAPACK `gesvd`）。
- 新增非原地版本 **`tsvd`**（矩阵版 `tsvd(a; trunc, alg)` 与张量版 `tsvd(a, left, right; ...)`），先复制输入再调用 `tsvd!`。
- 张量版 `tsvd!(a, left, right; ...)`：由于维度置换必须内部拷贝，输入张量本身**不会被修改**；`!` 表示复用内部拷贝作为 workspace（`leftorth!`/`rightorth!` 张量版同理，文档字符串已如实注明）。

### `leftorth` / `rightorth`（新增）

- 新增非原地版本并 export：`leftorth(A; alg, atol)`、`leftorth(A, left, right; ...)`、`rightorth(A; alg, atol)`、`rightorth(A, left, right; ...)`，复制输入后调用对应的 `!` 版本。
- `!` 版本的 `alg` 取值不变：`leftorth!` 支持 `QR`/`QRpos`/`SVD`/`SDD`/`Polar`，`rightorth!` 支持 `LQ`/`LQpos`/`SVD`/`SDD`/`Polar`。

### `isometry`（重命名并 export）

- `_eye` 重命名为 **`isometry`**，全部调用点已更新：
  - `isometry(T, m, n)` / `isometry(T, d)`（与原 `_eye` 相同）；
  - 新增便捷方法 `isometry(m, n)`、`isometry(d)`（默认 `Float64`）。

### 删除的未使用函数

- `texp`、`move_selected_index_forward`、`move_selected_index_backward`：无任何调用点，直接删除。
- `easy_swap!` / `naive_swap!` / `_swap_gate`：死代码（`@tensor` 的 `A[a,b;c,d]` 分号语法只产生普通 4D 数组，内部对 4D 张量调用仅接受矩阵的 `tsvd!` 本会 MethodError，且 `permute!(ADT, ...)` 从未被调用），已注释掉。

## export 变更

```julia
# 新增 export
tsvd, leftorth, rightorth, isometry
# 移除 export
stable_tsvd!
```

## 测试

- `test/auxiliary.jl` 新增 `isometry`（8 项）、`tsvd`（23 项，覆盖 SVD/SDD 驱动、输入不被修改、截断、张量版本）、`leftorth and rightorth`（47 项，覆盖全部算法、atol 截断、张量版本）三个 testset；原依赖"`tsvd!` 不修改输入"的用例改用 `tsvd`。
- 注意：MatrixAlgebraKit 的 `rightorth!(A, Polar(), ...)` 要求宽矩阵（列数 ≥ 行数），tall 矩阵会抛 `ArgumentError`。
- 全套测试通过（93 个 testset，0 Fail / 0 Error）。
