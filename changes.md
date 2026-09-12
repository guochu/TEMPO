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
