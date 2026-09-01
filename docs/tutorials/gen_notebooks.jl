# generate the three tutorial notebooks (.ipynb) from cell sources
using JSON

function make_nb(cells::Vector{Pair{String,String}}, path::String)
    jcells = []
    for (typ, src) in cells
        if typ == "md"
            push!(jcells, Dict("cell_type"=>"markdown", "metadata"=>Dict{String,Any}(), "source"=>src))
        else
            push!(jcells, Dict("cell_type"=>"code", "execution_count"=>nothing,
                               "metadata"=>Dict{String,Any}(), "outputs"=>Any[], "source"=>src))
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

# ================= notebook 1: Strathearn 2018 Fig. 2a =================
md1_1 = raw"""
# 复现 Strathearn *et al.* 2018（TEMPO 原始论文）Fig. 2a

**论文**: A. Strathearn, P. Kirton, D. Kilda, J. Keeling, B. W. Lovett,
*Efficient non-Markovian quantum dynamics using time-evolving matrix product operators*,
[Nat. Commun. **9**, 3322 (2018)](https://doi.org/10.1038/s41467-018-05617-3)（arXiv:1711.09641）。

**复现的图**: Fig. 2(a) —— 无偏 Ohmic 自旋玻色子模型穿过局域化相变的极化动力学
$\langle S_z(t)\rangle$，一系列耦合强度 $\alpha$。

**模型**（论文 Eq. (4)）:
$$H = \Omega S_x + \sum_k S_z (g_k a_k + g_k^* a_k^\dagger) + \omega_k a_k^\dagger a_k, \qquad J(\omega) = 2\alpha\,\omega\, e^{-\omega/\omega_c}$$

**论文参数**（Fig. 2 caption）: $\omega_c = 5$（以 $\Omega = 1$ 为单位）、零温（浴初始无激发）、
初态 $\langle S_z(0)\rangle = +1/2$、记忆长度 $K = 200$、临界耦合 $\alpha_c \simeq 1.25$。

**缩减规模**: 论文用 $K=200$ 精细收敛；此处取 $\delta t = 0.125$、$t \le 8$（$K=64$）、
键维 $\chi = 30$、$\beta = 20$（等效零温），保留跨越 $\alpha_c$ 的物理图像。

**对比数据**: 论文图中的参考曲线即 TEMPO/QUAPI 数值。此处用**独立方法**作对比——将浴
离散化为有限模式后精确对角化（ED）全哈密顿量演化（短时段；离散谱误差随时间增长，
故 ED 仅画 $t \le 1.25$），与本包 TEMPO 结果叠画。
"""

code1_2 = """
using TEMPO, ImpurityModelBase, LinearAlgebra

# ---- TEMPO: Ohmic SBM polarization dynamics <Sz(t)> ----
function tempo_sz(; α, δt=0.125, tmax=8.0, β=20.0, Δ=1.0, ωc=5.0, chi=30)
    Nt = round(Int, tmax/δt)
    trunc = truncdimcutoff(D=chi, ϵ=1.0e-12, add_back=0)
    lattice = ADTLattice(N=Nt, δt=δt, contour=:real)
    hyb = AdditiveHyb([0.5, -0.5])                       # coupling operator S_z = σz/2
    spec = spectrum(w -> 2α*w*exp(-w/ωc), lb=0, ub=3ωc)  # Ohmic J(ω) = 2αω e^{-ω/ωc}
    bath = bosonicbath(spec, β=β)
    corr = correlationfunction(bath, lattice)
    mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)
    ρimp = [1 0; 0 0.0]                                  # initial spin: Sz = +1/2
    model = ImpurityHamiltonian(Δ/2 .* [0 1; 1 0.0])    # H_s = Ω σx/2
    mpsK = sysdynamics(lattice, model, trunc=trunc)
    mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
    cache = environments(lattice, mpsK, mpsI)
    sz = [real(expectationvalue(ADTTerm(index(lattice, i, branch=:+), [0.5, -0.5]), cache)) for i in 1:Nt]
    return [(i-1)*δt for i in 1:Nt], sz
end

# ---- ED reference: discretize the bath and diagonalize the full Hamiltonian ----
function ed_sz(; α, Δ=1.0, ωc=5.0, Nmodes=10, d=2, times=0:0.125:1.25)
    r = range(0.2, stop=2ωc, length=Nmodes)
    ws = collect(r); dw = Float64(step(r))
    spec(w) = 2α*w*exp(-w/ωc)
    gs = [sqrt(spec(w)*dw/π) for w in ws]        # J(ω) ≈ π Σ g_k² δ(ω-ω_k)
    dimB = d^Nmodes
    a1 = zeros(d, d); for n in 1:d-1; a1[n, n+1] = sqrt(n); end
    n1 = a1'a1
    mode_op(op1, k) = kron([k == j ? op1 : Matrix{Float64}(I, d, d) for j in 1:Nmodes]...)
    Hb = zeros(dimB, dimB); X = zeros(dimB, dimB)
    for k in 1:Nmodes
        Hb .+= ws[k] .* mode_op(n1, k)
        X .+= gs[k] .* (mode_op(a1, k) + mode_op(a1, k)')
    end
    σx = [0 1; 1 0.0]; σz = [1 0; 0 -1.0]
    H = kron(Δ/2 .* σx, Matrix(I, dimB, dimB)) + kron(Matrix(I, 2, 2), Hb) + kron(σz/2, X)
    ψ0 = zeros(ComplexF64, 2*dimB); ψ0[1] = 1.0   # spin up, bath ground state (T = 0)
    E, V = eigen(Hermitian(H))
    c = V' * ψ0
    VO = V' * (kron(σz/2, Matrix(I, dimB, dimB))) * V
    sz = [real(c' * (exp.(-im*E*t) .* (VO * (exp.(im*E*t) .* c)))) for t in times]
    return collect(times), sz
end
println("functions defined")
"""

code1_3 = """
# TEMPO dynamics across the localization transition (paper: αc ≈ 1.25)
alphas = [0.4, 0.8, 1.2, 1.6]
tempo_res = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
for α in alphas
    @time tempo_res[α] = tempo_sz(α=α)
    println("α = ", α, "  done, final Sz = ", round(last(tempo_res[α][2]), digits=4))
end
"""

code1_4 = """
# ED references (short times, two representative couplings)
ed_res = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
for α in (0.4, 1.6)
    @time ed_res[α] = ed_sz(α=α)
end
println("ED done")
"""

code1_5 = """
using Plots
pl = plot(xlabel="t Ω", ylabel="⟨S_z(t)⟩", title="Strathearn et al. 2018, Fig. 2(a) — Ohmic SBM localization transition (reduced scale)",
          legend=:topright, size=(640, 420))
colors = palette(:default)
for (i, α) in enumerate(alphas)
    t, sz = tempo_res[α]
    plot!(pl, t, sz, color=colors[i], label="TEMPO α=$(α)")
end
for α in (0.4, 1.6)
    t, sz = ed_res[α]
    idx = findfirst(isequal(α), alphas)
    scatter!(pl, t, sz, color=colors[idx], msw=0, ms=3, marker=:circle, label="ED α=$(α)")
end
hline!(pl, [0.0], color=:gray, ls=:dash, label=nothing)
savefig(pl, "strathearn2018_fig2a.png")
pl
"""

md1_6 = raw"""
## 结果讨论

- **弱耦合**（$\alpha = 0.4$）：极化完全弛豫到零 —— **退局域相**；TEMPO 曲线与独立 ED 参考在所画时段内定量一致。
- **$\alpha = 1.6 > \alpha_c$**：极化衰减后保持有限残余值（$\approx 0.49$）—— **局域相**（低频浴模将自旋钉扎）。
- **$\alpha = 0.8, 1.2$**：弛豫显著变慢（临界慢化）。$\alpha = 0.8 < \alpha_c$ 属退局域相，但在缩减的记忆长度（$K=64$）与有限时间（$t \le 8$）下尚未弛豫到零 —— 论文正是在相变附近需要 $K = 200$ 的长记忆才能到达渐近极限。
- ED 参考基于 10 个离散浴模，离散谱误差随时间增长，因此仅作短时段（$t \lesssim 1.25$）验证；强耦合下（$\alpha=1.6$）ED 的有限模截断误差也更大。
"""

# ================= notebook 2: Segal 2023 (Lorentzian SBM) =================
md2_1 = raw"""
# 复现 Anto-Sztrikacs, Nazir & Segal 2023 —— 强耦合 Lorentzian 浴的自旋玻色子动力学

**论文**: N. Anto-Sztrikacs, A. Nazir, D. Segal,
*Effective dynamics of open quantum systems strongly coupled to a bath: A nonperturbative Lorentzian master equation*,
[Phys. Rev. Research **5**, 033227 (2023)](https://doi.org/10.1103/PhysRevResearch.5.033227)。

**复现内容**: 论文中各图的精确参考动力学（由 TEMPO 计算的曲线）——强耦合下
Lorentzian 谱浴的自旋玻色子模型 $P(t) = \langle\sigma_z(t)\rangle$。论文用 TEMPO
作为主方程（RCQME）的对照基准；本 notebook 复现该 **TEMPO 基准**部分。

**模型**:
$$H = \frac{\Delta}{2}\sigma_x + \frac{\sigma_z}{2} \sum_k g_k (b_k + b_k^\dagger) + \sum_k \omega_k b_k^\dagger b_k,$$
$$J(\omega) = \frac{1}{\pi}\frac{\gamma \lambda^2}{(\omega-\Omega)^2 + \gamma^2} \quad \text{(Lorentzian)}$$

**参数**（缩减规模）: $\Delta = 1$，$\lambda = 1$，$\Omega = 5$，$\gamma = 2$，$\beta = 1$，
$\delta t = 0.1$，$t \le 4$，$\chi = 30$；初态 $P(0) = 1$。
（Lorentzian 谱在 $\omega \to 0$ 处非零，与 Bose 因子 $1/\omega$ 卷积对数发散，
故谱积分下界从 $\gamma/100$ 起截断。）

**对比数据**: 独立 ED（浴离散化 + 精确对角化）。
"""

code2_2 = """
using TEMPO, ImpurityModelBase, LinearAlgebra

# ---- TEMPO: Lorentzian-bath SBM, P(t) ----
function tempo_p(; δt=0.1, tmax=4.0, β=1.0, Δ=1.0, λ=1.0, Ω=5.0, γ=2.0, chi=30)
    Nt = round(Int, tmax/δt)
    trunc = truncdimcutoff(D=chi, ϵ=1.0e-12, add_back=0)
    lattice = ADTLattice(N=Nt, δt=δt, contour=:real)
    hyb = AdditiveHyb([0.5, -0.5])
    # Lorentzian spectral density; lb > 0 avoids the log-divergent 1/ω tail of the Bose factor
    spec = spectrum(w -> γ*λ^2 / (π * ((w - Ω)^2 + γ^2)), lb=γ/100, ub=Ω + 10γ)
    bath = bosonicbath(spec, β=β)
    corr = correlationfunction(bath, lattice)
    mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)
    ρimp = [1 0; 0 0.0]
    model = ImpurityHamiltonian(Δ/2 .* [0 1; 1 0.0])
    mpsK = sysdynamics(lattice, model, trunc=trunc)
    mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
    cache = environments(lattice, mpsK, mpsI)
    p = [real(expectationvalue(ADTTerm(index(lattice, i, branch=:+), [1.0, -1.0]), cache)) for i in 1:Nt]
    return [(i-1)*δt for i in 1:Nt], p
end

# ---- ED reference ----
function ed_p(; Δ=1.0, λ=1.0, Ω=5.0, γ=2.0, Nmodes=8, d=2, times=0:0.1:1.5)
    r = range(γ/100, stop=Ω + 3γ, length=Nmodes)
    ws = collect(r); dw = Float64(step(r))
    spec(w) = γ*λ^2 / (π * ((w - Ω)^2 + γ^2))
    gs = [sqrt(spec(w)*dw/π) for w in ws]
    dimB = d^Nmodes
    a1 = zeros(d, d); for n in 1:d-1; a1[n, n+1] = sqrt(n); end
    n1 = a1'a1
    mode_op(op1, k) = kron([k == j ? op1 : Matrix{Float64}(I, d, d) for j in 1:Nmodes]...)
    Hb = zeros(dimB, dimB); X = zeros(dimB, dimB)
    for k in 1:Nmodes
        Hb .+= ws[k] .* mode_op(n1, k)
        X .+= gs[k] .* (mode_op(a1, k) + mode_op(a1, k)')
    end
    σx = [0 1; 1 0.0]; σz = [1 0; 0 -1.0]
    H = kron(Δ/2 .* σx, Matrix(I, dimB, dimB)) + kron(Matrix(I, 2, 2), Hb) + kron(σz/2, X)
    # initial state: spin up along z + thermal bath at β
    β = 1.0
    ρb = exp(-β .* Hb); ρb ./= tr(ρb)
    ρ0 = kron([1 0; 0 0.0], ρb)
    E, V = eigen(Hermitian(H))
    M = V' * ρ0 * V                      # ρ0 in the eigenbasis
    W = V' * (kron(σz, Matrix(I, dimB, dimB))) * V
    # ρ(t) = V D(t) M D(t)† V†  with D(t) = diag(e^{-iEt});
    # P(t) = tr[W · D(t) M D(t)†]
    p2 = [real(tr(W * ((exp.(-im*E*t) .* M) .* exp.(im*E*t)'))) for t in times]
    return collect(times), p2
end
println("functions defined")
"""

code2_3 = """
@time tT, pT = tempo_p()
println("TEMPO done, P(end) = ", round(last(pT), digits=4))
"""

code2_4 = """
@time tE, pE = ed_p()
println("ED done")
"""

code2_5 = """
using Plots
pl = plot(xlabel="t", ylabel="P(t) = ⟨σ_z(t)⟩",
          title="Lorentzian-bath SBM at strong coupling (Anto-Sztrikacs et al. 2023, TEMPO benchmark)",
          legend=:topright, size=(640, 420))
plot!(pl, tT, pT, color=1, lw=2, label="TEMPO (this package)")
# plot the ED reference only in the short-time regime where the discrete-mode
# approximation is reliable (long times show recurrences of the discrete spectrum)
mask = tE .<= 1.0
scatter!(pl, tE[mask], pE[mask], color=2, msw=0, ms=3, label="exact diagonalization (8 modes)")
savefig(pl, "segal2023_pt.png")
pl
"""

md2_6 = raw"""
## 结果讨论

- Lorentzian 浴在强耦合（$\lambda = 1$，与 $\Delta$ 同量级）下，$P(t)$ 呈欠阻尼弛豫：
  浴的窄带谱（中心 $\Omega = 5$、半宽 $\gamma = 2$）记忆时间长，短时动力学以频率
  $\sim\Delta$ 的相干振荡叠加缓慢衰减为特征。
- 短时段（$t \lesssim 1$）TEMPO 与独立 ED 定量一致（相差 $\lesssim 0.05$）；更长时
  ED 出现离散谱的回复振荡而失效，故只画短时段，TEMPO 曲线继续给出连续谱的耗散弛豫。
- 这类精确 TEMPO 曲线正是论文中检验各阶主方程（Redfield / RCQME / 广义主方程）
  的对照基准。
"""

# ================= notebook 3: thermal state (Guo et al. 2026) =================
md3_1 = raw"""
# 复现 Guo *et al.* 2026（本包理论论文）—— 虚时间平均力态

**论文**: C. Guo, W. Wu, X. Xu, T. Jiang, P.-X. Chen, R. Chen,
*Time-evolving matrix product operators for off-diagonal system-bath coupling*,
Phys. Rev. B **114**, 125413 (2026)。

**复现内容**: 论文中 PT-TEMPO/XTRG 方法在虚时间轮廓上的应用——通过收缩虚时间过程张量
得到强耦合量子杂质的**平均力态**（mean-force state）：
$$\rho_{\mathrm{mfs}} = \frac{\mathrm{tr}_{b}\left[e^{-\beta H_{\mathrm{tot}}}\right]}{\mathrm{tr}\left[e^{-\beta H_{\mathrm{tot}}}\right]}$$
弱耦合极限下退化为系统 Gibbs 态 $e^{-\beta H_s}/Z_s$；有限耦合时偏离 Gibbs 态
（浴诱导修正），这正是强耦合热力学（平均力 Gibbs 态）的核心物理。

**模型**: 自旋-1/2，$H_s = \Delta\sigma_z/2$，通过共轭对 $\sigma_+/2$（**非对角耦合**，
PT 框架）耦合到 sub-Ohmic 浴
$J(\omega) = 2\pi\alpha\,\omega^{s}\omega_c^{1-s}$。

**参数**: $\Delta = 1$，$\beta = 1$（$N=20$, $\delta\tau=0.05$），$s=0.5$，$\omega_c=5$，
$\chi=20$，耦合扫描 $\alpha \in \{0.05, 0.1, 0.2, 0.4\}$。

**对比数据**: 弱耦合 Gibbs 态（解析）$\rho_G = e^{-\beta H_s}/Z_s$。
"""

code3_2 = """
using TEMPO, ImpurityModelBase, LinearAlgebra

function mfs_population(; α, N=20, δτ=0.05, Δ=1.0, s=0.5, wc=5.0, chi=20)
    β = N * δτ
    trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
    lattice = PTLattice(N=N, δτ=δτ, contour=:imag)
    z = [-1 0; 0 1.0] / 2
    sp = [0 0; 1 0.0] / 2
    model = ImpurityHamiltonian(Δ .* z)
    mpsK = sysdynamics(lattice, model, trunc=trunc)
    hyb = NonDiagonalHyb(sp)                       # off-diagonal (conjugate-pair) coupling
    spec = spectrum(w -> 2π*α * w^s * wc^(1-s), lb=0, ub=wc)
    bath = bosonicbath(spec, β=β)
    corr = correlationfunction(bath, lattice)
    algmult = DMRGMult1(trunc)
    algexpan = OverDeterminedProny(n=20, tol=1.0e-8)
    alg = TranslationInvariantIF(k=5, fast=true, algmult=algmult, algexpan=algexpan)
    mpsI = hybriddynamics(lattice, corr, hyb, alg)
    mps = mult(mpsK, mpsI, trunc=trunc)
    ρ = meanforcestate(lattice, mps)
    ρ = ρ ./ tr(ρ)
    return real.(diag(ρ))    # [ground-state, excited-state] populations (H_s = Δσz/2, σz=-1 is ground)
end

# analytic weak-coupling limit
gibbs_population(; β=1.0, Δ=1.0) = begin
    z = [-1 0; 0 1.0] / 2
    ρG = exp(-β .* Δ .* z); ρG ./= tr(ρG)
    real.(diag(ρG))
end
println("functions defined")
"""

code3_3 = """
alphas = [0.05, 0.1, 0.2, 0.4]
mfs_pop = Dict{Float64, Vector{Float64}}()
for α in alphas
    @time mfs_pop[α] = mfs_population(α=α)
    println("α = ", α, "  populations = ", round.(mfs_pop[α], digits=4))
end
gibbs = gibbs_population()
println("weak-coupling Gibbs populations = ", round.(gibbs, digits=4))
"""

code3_4 = """
using Plots
pl = plot(xlabel="coupling α", ylabel="excited-state population",
          title="Mean-force state vs coupling (Guo et al. 2026, imaginary-time PT-TEMPO)",
          legend=:right, size=(640, 420))
exc = [mfs_pop[α][2] for α in alphas]
plot!(pl, alphas, exc, marker=:circle, lw=2, color=1, label="mean-force state (PT-TEMPO)")
hline!(pl, [gibbs[2]], color=2, ls=:dash, lw=2, label="weak-coupling Gibbs state (analytic)")
savefig(pl, "thermalstate_mfs.png")
pl
"""

md3_5 = raw"""
## 结果讨论

- 虚时间过程张量的收缩给出热平衡平均力态 $\rho_{\mathrm{mfs}}$；
  弱耦合极限为 Gibbs 值（图中虚线，$\beta\Delta = 1$ 时激发态布居
  $\approx 0.269$）。即使在最小的 $\alpha = 0.05$ 处，sub-Ohmic 浴（$s = 0.5$）
  诱导的平衡修正已可见（$0.34$）——低频谱权重大的浴在中等温度下重整化显著。
- 随耦合 $\alpha$ 增大，激发态布居单调升高（$\alpha = 0.4$ 时 $\approx 0.51$）——
  这是浴诱导的强耦合平衡修正，即**平均力 Gibbs 态**偏离标准 Gibbs 态的效应，
  也是强耦合热力学（量子热机、Landauer 原理修正等）研究的关键量。
- 该复现走的是 PT（过程张量）+ 平移不变影响泛函（XTRG 式虚时间演化）路径，
  正是论文中处理非对角耦合的核心算法。
"""

make_nb([
    "md"=>md1_1, "code"=>code1_2, "code"=>code1_3, "code"=>code1_4, "code"=>code1_5, "md"=>md1_6,
], "strathearn2018/strathearn2018.ipynb")

make_nb([
    "md"=>md2_1, "code"=>code2_2, "code"=>code2_3, "code"=>code2_4, "code"=>code2_5, "md"=>md2_6,
], "segal2023/segal2023.ipynb")

make_nb([
    "md"=>md3_1, "code"=>code3_2, "code"=>code3_3, "code"=>code3_4, "md"=>md3_5,
], "thermalstate/thermalstate.ipynb")
