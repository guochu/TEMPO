# TEMPO.jl

本工具包是 **时间演化矩阵乘积算符（Time-Evolving Matrix Product Operator, TEMPO）** 方法的一个 Julia 实现，其算法理论基于文献：

> C. Guo, W. Wu, X. Xu, T. Jiang, P.-X. Chen, and R. Chen,
> *Time-evolving matrix product operators for off-diagonal system-bath coupling*,
> **Phys. Rev. B 114, 125413 (2026)**.

与只支持对角系统-浴耦合的原始 TEMPO 不同，本实现基于**过程张量（Process Tensor, PT）**框架，将 TEMPO 推广到了更一般的**非对角系统-浴耦合**情形（系统通过一对共轭的非厄米算符 $\hat{A}^\dagger, \hat{A}$ 与浴耦合），并且统一了：

- 标准 TEMPO（对角、对易耦合，`ADT` + 部分影响泛函）；
- 对角但非对易的多浴耦合；
- 非对角（共轭对）耦合（`PT` + 平移不变影响泛函）；
- 实时间、虚时间、以及混合（Kadanoff–Baym）轮廓上的演化；
- 含时系统-浴耦合。

## 文档导航

| 页面 | 内容 |
|---|---|
| [快速上手](@ref) | 三类典型问题的完整计算流程（可直接照抄运行） |
| [使用手册](@ref) | 核心组件、超参数与误差来源、代码结构、文献对应 |
| [实践指南](@ref) | 实战经验：路径选择、收敛检验、常见坑与诊断（**做一般计算任务前建议先读**） |
| [实现细节](@ref) | 内部数据结构、算法流程与张量约定 |
| [API 参考](@ref) | 全部导出符号的逐项文档 |

配套的论文复现 notebook 见 `docs/tutorials/`（strathearn2018、otterpohl2025、guo2026、spinboson）。

## 安装

本工具包依赖以下包（见 `Project.toml`）：

| 包 | 作用 |
|---|---|
| `ImpurityModelBase` | 定义谱密度、浴（`bosonicbath`）、玻色算符等基础对象 |
| `QuAPI` | 提供 `ContourIndex`、`branch`、`index` 等轮廓基础类型 |
| `ExpExp` | 混合化函数的指数展开（Prony 方法） |
| `TensorOperations` | 张量网络缩并 |
| `Polynomials`、`LinearAlgebra`、`Statistics`、`TupleTools`、`Logging` | 通用数值与线性代数工具 |

在 Julia 中激活项目后即可使用：

```julia
using TEMPO, ImpurityModelBase, LinearAlgebra
```

!!! note
    `spectrum`、`bosonicbath`、`bosonaoperator`、`bosondensityoperator`、`Leggett` 等函数定义在 `ImpurityModelBase` 中，需要一并 `using`。

## 核心概念

### 量子杂质问题（QIP）

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

### Feynman–Vernon 影响泛函（IF）

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

### 两种张量网络：ADT 与 PT

| 对象 | 全称 | 表示 | 适用耦合 |
|---|---|---|---|
| `ADT` | 增强密度张量（Augmented Density Tensor） | MPS | 对角耦合（原始 TEMPO） |
| `ProcessTensor` (`PT`) | 过程张量 | MPO | 对角非对易、非对角耦合（文献扩展） |

一个 PT 可通过对相邻位点施加 3D copy 张量系统地转换为 ADT（文献 Fig. 4）。

### 轮廓（Contour）

本工具包支持三种时间轮廓，通过 `contour` 关键字选择：

- `contour=:real`（等价于 `:Keldysh`）：实时间演化，系统初态为 $\hat{\rho}_S \otimes \hat{\rho}_B$；
- `contour=:imag`：虚时间演化（$0\to\beta$），对应有限温度配分函数与 Matsubara 关联函数；
- `contour=:mixed`（等价于 `:Kadanoff`）：L 形 Kadanoff–Baym 轮廓（虚时间 + 实时间混合）。
