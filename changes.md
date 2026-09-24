# 跟进（2026-09-23）：复用 FiniteMPSAlgorithms 新增接口

FMA 接口调整（`SchurMPOTensor` 构造器只接受完整逻辑块矩阵（m, n ≥ 2，链边界由
`tompotensors` 下游处理；单类型参数）等）后，清理 TEMPO 侧剩余的重复实现。全量
测试通过（1244/1244）：

- **mult 的 SVD 压缩路线整体委托**：删除 `adt/mult/svdmult.jl` 与
  `pt/mult/svdmult.jl`（zip-up 的 QR 累积 + 截断右扫本地实现，约 90 行）。
  ADT（Hadamard 乘积）改用 FMA 的 `hadamard(x.parent, y.parent, SVDCompression)`、
  PT（MPO 乘积）改用 FMA 的 `mult(x.parent, y.parent, SVDCompression)`——均为
  "精确乘积 + 单次 SVD 压缩扫"，scaling 语义一致（外标度经输出 `scaling` 字段
  继承）。`mult` / `mult!` 的 `trunc`/`verbosity` keyword 接口与
  `SVDCompression`/`DMRGAlgorithm` 分派保持不变（SVD 路线的 `verbosity` 由 FMA
  忽略）。顺带删除失去调用者的 `_mult_site_n`。
- **`_rescaling!` 删除**（`adt/orth.jl`）：与 FMA 的 `_renormalize!` /
  `_rescaling!` 逐行等价，`iterativemult` 的调用点改为
  `_renormalize!(z, z[1], false)`。
- **`randomadt` / `randompt` 保留本地实现**：FMA 的 `randommps` / `randommpo`
  用 `max_bonddims`（按单物理维）截顶键维并规范化，会改变"每键 = D"的构造
  语义（测试与调用方依赖它构造键维充足的 ALS ansatz），故不委托。
- **`SchurMPOTensor` 适配确认**：TEMPO 的衰减项 cell 恒为 (N+2)×(N+2)（N ≥ 0），
  满足新接口的 m, n ≥ 2 要求；`MPOHamiltonian([mpoj, ...])` + `tompotensors`
  的用法与新边界处理方式兼容，无需改动。
- 测试：谱对比 testset 的 DMRG1 配置 `maxiter` 提至 20（个别随机实例 5 轮
  ALS 不收敛，避免收敛噪声干扰 finalize 正确性检验）。
- `docs/src/internals.md` §10.2/§10.3 同步更新（SVD 路线委托与 finalize 描述）。

# Bug 修复（2026-09-23）：`mult`（DMRG1 路线）finalize 的 Schmidt 谱不正确

## 问题

`_finalize!`（`fmabackend.jl`，ALS 收敛后的收尾 sweep）的旧实现是在**环境规范系**下
对 ALS 局部目标 `mpsj = L·W·R` 做逐站点 tsvd，把局部目标的谱写入 `z.s`，且 `u·s`
只在 site==2 折回（其余站点直接丢弃）。这样得到的 `z.s` 不是最终输出态各键的
Schmidt 值：它只在 ALS 精确不动点、双侧规范等距、**且无截断**时才与真谱一致
（此时实测偏差 ~1e-16）；一旦 finalize 发生截断（如 `truncdim(D)` 且 D 小于精确
键维），存谱与输出态的偏差可达截断误差量级（实测 ~1e-3），`iscanonical` 检验
（`Diagonal(s²) ≈` 左环境）失败。此前测试未检出是因为 mult 只检验 `distance`。

## 修复

- **`_finalize!` 重写**（ADT/PT 共用一个方法）：ALS 收敛后输出链本身就是最终
  乘积态，规范与谱不再依赖环境——直接调用 FMA 的 `canonicalize!`（QR 左扫精确、
  不改变状态；SVD 右扫按 `trunc` 截断、每步把 `u·s` 折回左侧并把谱写入该键）。
  两步扫后链右正交，`z.s` 是（截断后）输出态各键的**精确** Schmidt 值，不依赖
  ALS 是否到达不动点；这与 FMA 自身 `_svd_mult` 的收尾方式一致。
  不再使用的 FMA 内部原语（`_reduce_hadamard_site` / `_reduce_site` /
  `_env_updateright` / `_updateright` / `_contract_last`）从绑定清单移除。
- **幅值记账修正**（`adt` / `pt` 的 `iterativemult.jl`）：`canonicalize!` 把链的
  范数因子折叠进 `scaling(z)`，wrapper 原来对 `scaling` 的绝对值赋值会丢掉该
  因子（输出幅值偏差实测达 7.6 倍）。改为与 svdmult 同构的乘法模式
  `setscaling!(z, scaling(z) * scaling(x) * scaling(y))`。

## 验证与测试

- 数值验证（6 站随机态，d=2，输入 D=4）：无截断（D=64）时与精确乘积的
  `distance ≈ 2.7e-8`、逐键谱偏差 ~1e-16、`iscanonical = true`；有截断（D=6）
  时 `iscanonical = true`（旧实现为 false）、谱为截断态的精确谱。
- `test/api/mps.jl`：multiplications 全部算法组合的输出新增 `@test iscanonical`；
  新增 `mult output: canonical form and spectra` testset——DMRG1 各 initguess
  路径的输出须 `iscanonical` 且无截断时逐键谱与精确乘积一致（允许近零 padding
  方向）；强截断（D=6）下输出须 `iscanonical`（其变分误差与 finalize 无关，
  新旧实现同量级）。
- 全套测试通过。

---

# 跟进（2026-09-22）：适配 FiniteMPSAlgorithms 接口变更

FMA 更新（scaling-safe 的 MPO 算术、链级 `swap!`/`permute!`/`permute`、
`Base.sum(::CanonicalMPS)`、`copyphydims` 等）后，TEMPO 侧相应简化，全部测试
继续通过（1181/1181）：

- **MPO 算术（scaling-safe）**：PT 的 `Base.:*` 直接用 FMA 的 CanonicalMPO 乘积
  （结果自带 `scaling(x)·scaling(y)`，删除手动挂 scaling）；PT / ADT 的
  `Base.:+` / `-` 直接委托 FMA（scaling 折入数据的块对角直和），删除本地实现与
  `L == 1` 特例。
- **链级置换**：ADT / PT 的 `swap!` / `permute!` / `permute` 改为对 payload 的
  委托包装（FMA 新增 CanonicalMPO 版 Hastings 式 `swap!` 与链级 `permute!`，
  未初始化 Schmidt 谱时自动规范化），删除本地 `_swap_gate` / `_permute!`；TEMPO
  侧的默认截断（`DefaultITruncation` / `DefaultKTruncation`）由包装保留。
- **integrate**：`integrate(x)` 委托 FMA 的 `sum(::CanonicalMPS)`（所有振幅之和，
  含 `scaling^L`）。`integrate(x, y)` **保持原 lazy 实现**——它是"x ⊙ y 的振幅
  之和"（物理指标共享求和、无共轭），与 FMA `dot`（内积）语义不同，不采用。
- **ADT×ADT 的 mult（DMRG1 路线）**：仍走 FMA 的 `HadamardCache`（物理指标共享
  的逐点乘积，比 `copyphydims` + `MultCache` 便宜且含义相同）；初始猜测
  `:svd` 改用 FMA 的 `svdguess_hadamard`（替代本地流式 SVD 猜测，仅猜测方向
  右到左的差异，ALS 收敛后结果一致）。finalize（QR 左扫 + 截断 SVD 右扫 + 键谱
  写入）共享化到 `fmabackend.jl` 的 `_finalize!`（`HadamardCache` / `MultCache`
  两个方法；PT 侧以 `NoTruncation()` 调用，保持旧行为）。PT×PT 的 mult（DMRG1
  路线）的 `:svd` 初始猜测同样改用 FMA 的 `svdguess_mult`（MPO×MPO 版流式
  SVD），删除本地的 `_svd_guess`。
- **import 修复**：`permute!` 此前不在 import 清单，TEMPO wrapper 成了独立新
  函数、委托调用无法命中 FMA 方法（落进 `Base.permute!` 兜底）——已加入 import。
  同时清理不再使用的导入（`permutation2swaps` / `_fused_pair` /
  `_contract_first`），新增 `svdguess_hadamard`。

# 重构（2026-09-22）：MPS/MPO 运算后端切换为 FiniteMPSAlgorithms

TEMPO 的张量层与 MPS/MPO 算法不再自行维护，统一委托给本地包
**FiniteMPSAlgorithms**（`Pkg.develop` 依赖，`Project.toml` 已登记）；TEMPO 的
公共接口（`ADT` / `ProcessTensor` / `mult` / `canonicalize!` / `DMRG1` /
`MPOHamiltonian` / `timeevompo` …）全部保持不变，在其实现之上套轻量级 wrapper。
全套测试通过（api + adtmodels + ptmodels，0 Fail / 0 Error）。

## 后端对应关系

| TEMPO 模块 | 原实现 | 现后端（FiniteMPSAlgorithms） |
|---|---|---|
| `src/tensorops/`（已删除） | 自维护（截断方案、tsvd!/leftorth!、tie/permute/isometry…） | 同名导出（`FiniteMPSAlgorithms` 的 tensorops 即 TEMPO 版 vendor），TEMPO re-export |
| `src/algorithms.jl` | 自定义 `Orthogonalize` / `SVDCompression` 结构体 | 直接采用 FMA 的 `Orthogonalize` / `SVDCompression` 类型；`SVDCompression(trunc; verbosity)` positional 构造与 `similar` 由 TEMPO 补充 |
| `src/mpohamiltonian/`（精简） | `AbstractSparseMPOTensor` / `SparseMPOTensor` / `SchurMPOTensor` / `MPOHamiltonian` / `tompotensors` / `w1w2.jl`（WI/WII/ComplexStepper/timeevompo） | 全部改用 FMA 的对应类型与函数；TEMPO 保留长程衰减项（`ExponentialDecayTerm` / `GenericDecayTerm` / `PowerlawDecayTerm`、`expand_decayterm`，其块矩阵构造 `SchurMPOTensor(cell)` 与 FMA 的 Schur 形式直接兼容）与 `compat.jl`（`TimeEvoMPOAlgorithm` 别名）。`phydim` 对稀疏 MPO 张量的方法直接由 FMA 提供（`import` 后共用同一泛型函数） |
| `src/adt/`、`src/pt/` 的 orth/linalg/mult | 自维护 QR/SVD sweep、ALS 引擎 | `leftorth!` / `rightorth!` / `canonicalize!` 委托 FMA 在内层 `CanonicalMPS`/`CanonicalMPO` payload 上的 `_leftorth!` / `_rightorth!` / `_canonicalize!`；ALS 引擎改用 FMA 的 `HadamardCache`（ADT×ADT：物理指标共享的逐点乘积）与 `MultCache`（PT×PT：MPO×MPO 乘积）+ `iterative_compute!` |

关键机制：

- **存储直用 FMA 链**：`ADT` / `ProcessTensor` 的存储 payload（字段 `.parent`）直接内嵌
  FiniteMPSAlgorithms 的 `CanonicalMPS` / `CanonicalMPO`（见 `adt/def.jl` / `pt/def.jl`）。
  FMA 的算法直接作用在 payload 上并就地改写，无需转换与拷回。`.data` 经 `getproperty`
  委托到 payload 的站点张量 `Vector`（与旧版语义一致）；`.s` / `.scaling` 委托到 payload
  的 Schmidt 值 / scaling（`propertynames` 同步声明，保证 `hasproperty` 守卫的 FMA 原语
  正常工作），全部消费方（`TransferMatrix`、`tdvpif`、`swap!` 等）无需改动。
- **接口变更**：`increase_bond!` 删除，由 FMA 的 `changebond!` 代替（TEMPO 提供
  `changebond!(psi::ADT, D::Int)` / `changebond!(g::ProcessTensor, D::Int)` 包装，委托
  payload 上的 `changebond!`）。与旧版只增不减不同，`changebond!` 会把键型整备到
  `min(D, feasible)`（超出部分按前导索引切片、不足部分补零），MPS 版随后无截断重新规范化。
- **单一泛型函数**：`mult`、`canonicalize!`、`leftorth!`、`rightorth!`、`swap!`、
  `distance`、`scaling`、`setscaling!`、`svectors_uninitialized`、`unset_svectors!`
  等通过 `import` 扩展，TEMPO 方法与 FMA 方法共存于同一函数。
- **DMRG1 翻译器**：TEMPO 的 `DMRG1(trunc; maxiter, tol, initguess, verbosity, callback)`
  保留原字段；驱动 FMA 引擎时经 `_fmadmrg1` 映射为 FMA 的 `DMRG1(maxiter, tol, D, verbosity)`。
- **finalize wrapper**：ALS 收敛后 TEMPO 特有的 finalize（QR 左扫 + 截断 SVD 右扫并把
  归一化键谱写入 `z.s`）保留在 TEMPO 侧，作用于 FMA 的缓存对象。

## 行为差异（有意保留 / 随之变化）

- `ADT` / `ProcessTensor` 的内部存储改为内嵌 FMA 的 `CanonicalMPS` / `CanonicalMPO`
  payload：对 `.data` / `.s` / `.scaling` 的字段访问经 `getproperty` 委托保持不变，
  但**旧版 `Serialization.serialize` 存出的 `ADT`/`ProcessTensor` 文件反序列化后类型
  不再匹配**（一次性影响，教程中缓存的 `.mps` 需重新生成）。
- `mult` 的 ALS 收敛判据改为 FMA 的 `iterative_compute!`（相邻 sweep 末位损失的相对变化，
  分母为 `|prev|`；旧实现分母为 `max(cur, prev)`），收敛结果等价，迭代数可能相差 1。
- `mult(x, y, DMRG1)` 的初始猜测仍由 TEMPO 侧流式 SVD 提供（`truncdim(alg.trunc.D)` 上限）；
  ALS 阶段不再截断，末次 finalize sweep 恢复 `alg.trunc` 截断（ADT）与键谱写入（ADT/PT，
  与旧版一致）。
- 无参 `SVDCompression()` 的默认截断为 FMA 的 `truncdimcutoff(D=64, ϵ=1e-12, add_back=0)`
  （原 TEMPO 默认 `D=100`）；显式传 `trunc` 的用法不受影响。其余 `Defaults`（D=100、tolgauge、
  DefaultMultAlg 等）不变。
- `XTRGIF` 与 `influenceoperatorstepper(s)` 的 `algmult` 类型约束放宽为 `MPSAlgorithm`
  （原 `DMRGAlgorithm`）：旧版 `SVDCompression` 是 `DMRGAlgorithm` 子类，切换到 FMA 类型后
  二者为平级算法，行为不变、类型约束放宽。

## 上游（FiniteMPSAlgorithms）同步修改

- `timeevompo(h::MPOHamiltonian{<:SchurMPOTensor}, dt)` 的二参便捷方法不再以
  `alg::MPSAlgorithm=WII()` 兜底（与 `ComplexStepper` 方法存在派发歧义），改为
  `timeevompo(h, dt) = timeevompo(h, dt, WII())`。

## 环境

- `Project.toml`：新增 `FiniteMPSAlgorithms`（develop），`[extras]`/`test` target 补上
  `Random`；移除不再直接使用的 `MatrixAlgebraKit`（现经 FMA 间接依赖）。

---

# Bug 修复（2026-09-17）：`swap!` 未交换物理指标

`swap!` / `permute!`（ADT 与 PT）的 swap gate 存在实现错误：两站合并张量的 SVD 分组为
`(l, p1) | (p2, r)`（物理指标留在各自一侧），重建出的链上物理指标顺序不变——swap 实际上是
一次恒等的规范重分解，没有交换任何物理指标。由于此前的 `permute` 测试只验证了自洽性
（`permute ∘ permute⁻¹ == id`、`iscanonical`），恒等实现同样能满足，因而未被检出。

## 修复

- 两站合并张量改为**先交换物理指标**再分解：ADT 按 `(l, p2, p1, r)`、PT（共轭对布局
  `(aL, pout, aR, pin)`）按 `(aL, pout2, pin2, pout1, pin1, aR)` 合并；SVD 分组后左因子
  携带 `bond+1` 站点的物理指标、右因子携带 `bond` 站点的物理指标，新键谱写回 `x.s[bond+1]`。
- 数值验证：单次相邻交换后两站块等于原始块在交换物理顺序下的结果（且不同于原始顺序）；
  相同相邻交换做两次完全还原；全链收缩（`integrate`）在交换前后保持不变。
- `test/api/mps.jl` 新增 `permute! moves the physical labels` testset：用 D=1 的 one-hot
  product 态（每个站点携带可区分的物理标签）直接验证置换后标签满足 `new[k] == old[perm[k]]`；
  原 `permute` testset 中「置换后仍 iscanonical」的断言移除（交换会把键谱移到其它键上，
  全键 Vidal 形式不再保持，物理正确性由「双交换还原 + 标签重排」两个断言保证）。

## 关联修复：`x.s` 的未初始化槽

排查中发现 `ADT{T,R}(data, scaling)` 与 `ProcessTensor{T,R}(data, scaling)` 内层构造器只初始化了
`s[1]` 与 `s[end]`，中间键谱槽是 `Vector{Union{Missing,Vector{R}}}` 的 **`undef` 引用**（访问即抛
`UndefRefError`），与「未定谱槽应为 `missing`」的设计（`unset_svectors!`、`svectors_uninitialized`
的 `any(ismissing, ...)`）不符。已修复：构造器中将 `s[2:end-1]` 显式初始化为 `missing`；
`svectors_uninitialized` 行为不变。

---

# 接口调整说明（2026-09-17）：`TruncateCutoff` / `trunccutoff` 更名

类型与构造函数更名，语义不变（按相对截断误差 ϵ 截断奇异值）。全套测试通过（0 Fail / 0 Error）。

| 旧名（已删除） | 新名 | 说明 |
|---|---|---|
| `TruncateCutoff` | `TruncateRelError` | 截断方案类型（含关键字构造 `TruncateRelError(; ϵ)`） |
| `trunccutoff` | `truncrelerr` | 便利构造函数（位置 `truncrelerr(ϵ)` 与关键字 `truncrelerr(; ϵ)` 等价） |

- 更名原因：与 `TruncateDim` / `TruncateDimCutoff` 命名对齐，明确 ϵ 的含义是**相对截断误差**而非绝对阈值。
- 迁移：`trunccutoff(ϵ)` → `truncrelerr(ϵ)`；类型注解 `::TruncateCutoff` → `::TruncateRelError`。
- `DefaultKTruncation = truncrelerr(Defaults.tolgauge)` 不变（仅名称）。

---

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
| `DefaultTruncation` | **已删除**（与 `DefaultITruncation` 合并；`mult!`/`canonicalize!`/`swap!` 的默认截断改用 `DefaultITruncation`） |
| `DefaultITruncation` | `truncdimcutoff(D=200, ϵ=1e-10)` → **`truncdimcutoff(D=Defaults.D, ϵ=Defaults.tolgauge)`**；`DefaultMultAlg = DMRG1(DefaultITruncation)` 随之变化 |
| `DefaultMultAlg` | 定义不变，取值跟随 `DefaultITruncation` |

替换位置：`boundarycondition!`、`_permute!`（ADT/PT linalg）、TTIIF 的 `_fit_to_lattice_diag/_offdiag`（adt/pt real）、TDVPIF 的 H 压缩（`_tdvpif_hamiltonian`）；`mult!`（zip-up SVD 乘法，ADT/PT svdmult.jl）、`canonicalize!`（ADT/PT orth.jl）、`swap!`（adt/def.jl、pt/linalg.jl）的默认截断由 `DefaultTruncation` 改为 `DefaultITruncation`。

## 其他删除

- `src/observables/correlations.jl`：整个文件删除（从未被 include 的死代码，`correlation` 函数无处分发使用；两点关联测量请用 `ADTTerm` 多点形式 + `apply!`/`integrate`，见 manual §Observables）。`docs/src/api.md`、`docs/src/manual.md` 同步清理。

## 迁移指南

- `SVDCompression(D=χ)` → `SVDCompression(truncdimcutoff(D=χ, ϵ=Defaults.tol))`（或按需写 `truncdimcutoff(D=χ, ϵ=...)`）。
- `SVDCompression(D=χ, tol=ε)` → `SVDCompression(truncdimcutoff(D=χ, ϵ=ε))`。
- `DMRG1(trunc)` / `DMRG1(trunc=...)`：`trunc` 需为 `truncdim(D)` 或 `truncdimcutoff(D, ϵ)`。
- `alg.D` / `alg.ϵ` 访问改为 `alg.trunc.D` / `alg.trunc.ϵ`。
- `DefaultIntegrationTruncation` / `DefaultMPOTruncation` / `DefaultTruncation` → `DefaultKTruncation` / `DefaultITruncation`。
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
