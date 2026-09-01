# generate otterpohl2025 notebook (.ipynb) from cell sources
using JSON

function make_nb(cells::Vector{Pair{String,String}}, path::String)
    jcells = []
    for (i, (typ, src)) in enumerate(cells)
        if typ == "md"
            push!(jcells, Dict("cell_type"=>"markdown", "metadata"=>Dict{String,Any}(), "id"=>"md$i", "source"=>src))
        else
            push!(jcells, Dict("cell_type"=>"code", "execution_count"=>nothing,
                               "metadata"=>Dict{String,Any}(), "id"=>"c$i", "outputs"=>Any[], "source"=>src))
        end
    end
    nb = Dict("cells"=>jcells,
              "metadata"=>Dict("kernelspec"=>Dict("display_name"=>"Julia 1.10", "language"=>"julia", "name"=>"julia-1.10"),
                               "language_info"=>Dict("name"=>"julia", "version"=>"1.10.11")),
              "nbformat"=>4, "nbformat_minor"=>5)
    open(path, "w") do f
        JSON.print(f, nb, 1)
    end
    println("wrote ", path)
end

md1 = raw"""
# 复现 Otterpohl *et al.* 2025：量子 $1/f^\eta$ 噪声驱动的自旋玻色弛豫

**论文**: F. Otterpohl, P. Nalbach, E. Paladino, G. A. Falci, M. Thorwart,
*Quantum 1/f$^\eta$ Noise Induced Relaxation in the Spin-Boson Model*,
[arXiv:2507.14329 (2025)](https://arxiv.org/abs/2507.14329)。

**复现的图**: Fig. 1(c) —— 谱指数 $s = -0.5$（即 $1/f^{1.5}$ 噪声区间）的自旋玻色模型在
$\alpha = 0.01,\,0.02,\,0.04$ 下的极化 $P(t) = \langle\sigma_z(t)\rangle$（实线）与
相干 $C(t) = \langle\sigma_x(t)\rangle$（虚线）。物理上：$s = -0.5$ 时
$\alpha = 0.01$ 仍处于相干区（衰减振荡），$\alpha = 0.02,\,0.04$ 进入**赝相干区**
（$P(t)$ 只剩单个极小，其位置由浴紫外截断 $1/\omega_c$ 的响应时间控制）——
与 Ohmic 情形（教程 strathearn2018）完全不同的弛豫机制。

**模型**（论文 Eq. (1)–(2)，$\hbar = k_B = 1$）:
$$H = \frac{\Omega}{2}\sigma_x + \frac{1}{2}\sigma_z \hat\xi, \qquad
J(\omega) = 2\alpha\, \omega_c^{1-s}\, \omega^{s}\, e^{-\omega/\omega_c}\,\Theta(\omega - \omega_{ir})$$
谱指数 $s=-1$ 对应 $1/f$ 噪声；本文取 $s = -0.5$（$1/f^{0.5}$），低频下 $J(\omega)\propto \omega^{-0.5}$
发散（重整化能发散），但浴关联函数 $C(t) \propto (1/\omega_c + it)^{-(s+1)}$ 仍收敛，
TEMPO 只需要关联函数，因此可以处理。

**论文参数**: $\Omega = 1$（能量单位）、$\omega_c = 10\Omega$、$\omega_{ir} = 0$、
$T = 0$（$\beta = \infty$）、初态因子化（$\langle\sigma_z(0)\rangle = +1$，浴为真空）、
时间范围 $\Omega t \le 5.14$。

**论文数据**: 从 arXiv PDF 矢量图数字化（`tutorial2_fig_data.json`），按 legend 颜色映射
（蓝 $\alpha=0.01$、橙 $0.02$、绿 $0.04$；实线 $P$、虚线 $C$）。

**本复现**: $\delta t = 0.05$、键维 $\chi = 40$、$\beta = \infty$、$t \le 5.15$。
$s<0$ 的时间幂律记忆使影响泛函键维随时间持续增长，$t\le 5$ 内 $\chi=40$ 已足够
（收敛性检验见结果讨论）。
"""

code2 = raw"""
using TEMPO, ImpurityModelBase, LinearAlgebra
using JSON, Serialization

# ---- TEMPO (PT path): 1/f^eta SBM, P(t) = <sigma_z>, C(t) = <sigma_x> ----
# NonDiagonalHyb(op) means op*a + op'*a' (op couples to the bath annihilation
# operator); for a real diagonal op this is identical to the additive coupling V(a+a').
# The PT lattice route is required so that the off-diagonal <sigma_x> is measurable.
function tempo_pc(; α, s=-0.5, δt=0.05, tmax=5.15, wc=10.0, chi=40, k=7, n=20)
    Nt = round(Int, tmax/δt)
    trunc = truncdimcutoff(D=chi, ϵ=1.0e-12, add_back=0)
    lattice = PTLattice(N=Nt, δt=δt, contour=:real)
    z = [1.0 0.0; 0.0 -1.0]                                          # σ_z (ContourOperator needs a matrix)
    x = [0.0 1.0; 1.0 0.0]                                           # σ_x
    hyb = NonDiagonalHyb(Matrix(Diagonal(0.5 .* [1.0, -1.0])))       # σ_z/2
    spec = spectrum(w -> 2α * wc^(1-s) * w^s * exp(-w/wc), lb=0, ub=5wc)
    bath = bosonicbath(spec, β=Inf)                                  # T = 0
    corr = correlationfunction(bath, lattice)
    mkpath("data")
    mpspath = "data/otterpohl2025_pt_betaInf_dt$(δt)_alpha$(α)_s$(s)_wc$(wc)_chi$(chi)_N$(Nt).mps"
    if ispath(mpspath)
        mpsI = Serialization.deserialize(mpspath)
    else
        algmult = DMRGMult1(trunc, maxiter=10)
        algexpan = OverDeterminedProny(n=n, tol=1.0e-8)
        alg = TranslationInvariantIF(k=k, fast=true, algmult=algmult, algexpan=algexpan)
        mpsI = hybriddynamics(lattice, corr, hyb, alg)
        Serialization.serialize(mpspath, mpsI)
    end
    ρimp = [1.0 0; 0 0.0]                                            # <σ_z(0)> = +1
    model = ImpurityHamiltonian(0.5 .* [0 1; 1 0.0])                 # H_s = Ω σ_x/2
    mpsK = sysdynamics(lattice, model, trunc=trunc)
    mps = mult!(mpsK, mpsI, trunc=trunc)
    cache = environments(lattice, mps, ρ₀=ρimp)
    p = [real(expectationvalue(ContourOperator(ContourIndex(i, :+), z), cache)) for i in 1:Nt]
    c = [real(expectationvalue(ContourOperator(ContourIndex(i, :+), x), cache)) for i in 1:Nt]
    return [(i-1)*δt for i in 1:Nt], p, c
end
"""

code3 = raw"""
paper2 = JSON.parsefile("tutorial2_fig_data.json")

tempo2_res = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}()
mkpath("result")
for α in (0.01, 0.02, 0.04)
    # chi=40 suffices for alpha=0.01/0.02 (checked: dt0.025/chi80 changes ~0.01);
    # alpha=0.04 requires chi=100 (chi=40 non-physical, chi=60 still unconverged)
    chi = (α == 0.04) ? 100 : 40
    @time tempo2_res[α] = tempo_pc(α=α, chi=chi)
    t, p, c = tempo2_res[α]
    println("α = ", α, "  done,  P(end) = ", round(p[end], digits=4), "  Cmin = ", round(minimum(c), digits=4))
    open("result/otterpohl2025_fig1c_alpha$(α)_s-0.5_wc10_dt0.05_chi$(chi).json", "w") do f
        write(f, JSON.json(Dict("t"=>t, "P"=>p, "C"=>c)))
    end
end
"""

code4 = raw"""
using Plots
pl = plot(layout=(1, 2), size=(960, 400),
          xlabel=["Ωt" "Ωt"], ylabel=["P(t)=⟨σz⟩" "C(t)=⟨σx⟩"],
          title=["s = -0.5: P(t)" "s = -0.5: C(t)"], legend=:bottomright)
for (i, α) in enumerate((0.01, 0.02, 0.04))
    # paper P (solid) and C (dashed), digitized
    ptsP = paper2["P_alpha$(α)_solid"]; ptsC = paper2["C_alpha$(α)_dashed"]
    plot!(pl[1], first.(ptsP), last.(ptsP), color=i, lw=2, label="paper α=$(α)")
    plot!(pl[2], first.(ptsC), last.(ptsC), color=i, lw=2, label="paper α=$(α)")
    # this package
    t, p, c = tempo2_res[α]
    plot!(pl[1], t, p, color=i, ls=:dash, lw=1.5, label="TEMPO α=$(α)")
    plot!(pl[2], t, c, color=i, ls=:dash, lw=1.5, label="TEMPO α=$(α)")
end
savefig(pl, "otterpohl2025_fig1c_vs_paper.png")
pl
"""

md5 = raw"""
## 结果讨论

图中实线为**从论文 Fig. 1(c) 数字化提取的原始数据**，虚线为本包 TEMPO 计算
（$\delta t = 0.05$、$\beta = \infty$、$t \le 5.15$；键维 $\chi = 40$
（$\alpha = 0.01/0.02$）与 $\chi = 100$（$\alpha = 0.04$））：

- **$P(t)$**：三条曲线都从 1 出发单调弛豫（$s=-0.5$ 时无相干振荡）
  - $\alpha = 0.01$：弛豫最快，$P$ 在 $t \approx 5$ 处仍在下降（相干区边界）；
  - $\alpha = 0.02, 0.04$：出现单个极小后缓慢回升（赝相干区特征），回升尺度由
    $1/\omega_c$ 浴响应时间控制，$\alpha$ 越大平台越高（浴重整化越强）。
- **$C(t)$**：全部变负且无过零振荡，$|C_{min}|$ 随 $\alpha$ 单调减小。
- 论文正文由此判定：$s=-0.5$ 时 $\alpha = 0.01$ 属相干区、$\alpha = 0.02/0.04$
  属赝相干区，与 Fig. 2 的动力学相图一致——本复现的曲线形状支持同样结论。

收敛性检验：$\alpha = 0.02$ 处将 $(\delta t, \chi) = (0.05, 40) \to (0.025, 80)$，
$P(t)$、$C(t)$ 的最大变化在 $\sim 0.01$ 量级。$\alpha = 0.04$ 处收敛要求显著更高：
$\chi = 40$ 的结果非物理（$P$ 超过 1、幅度漂移），$\chi = 60$ 仍未收敛（与
$\chi = 100$ 差异可达 $\sim 0.5$），最终采用 $\chi = 100$。$s < 0$ 使浴关联函数呈
幂律长尾（$C(t) \propto t^{-0.5}$），影响泛函键维随 $t$ 持续增长且耦合越强增长
越快——这是本算例与 Ohmic 教程（指数衰减记忆）的本质区别。
"""

make_nb([
    "md" => md1,
    "code" => code2,
    "code" => code3,
    "code" => code4,
    "md" => md5,
], joinpath(@__DIR__, "otterpohl2025.ipynb"))
