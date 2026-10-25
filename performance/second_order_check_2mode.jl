# =============================================================================
# 二阶 QuAPI/TEMPO 离散化验证 —— 双模玻色浴（排除单模巧合的可能）
#
# 数学预期: 多模高斯浴逐模分解, η 核 = 各模核之和; IF 仍是分段常值路径的精确
# Feynman-Vernon 偏迹, 故端点半步旋转的二阶结论不应依赖模数。
# 本脚本用两个非对称、非公度的模 (ω₁,g₁)=(1.3,0.35), (ω₂,g₂)=(0.7,0.22) 检验:
#   (A) IF(MPS) vs 双模解析 η (N=3)
#   (B) 固定 T=1.6, 逐次减半 δt: 一阶 ~ δt¹, 二阶 ~ δt² (vs 双模全空间 ED)
#   (C) 二阶结果 ≈ 全空间 Strang, 一阶结果 ≈ 全空间 Lie (恒等式)
#
# 双模浴通过原生 DiscreteSpectrum 谱入口进入包（bosonic_Δt 经 DeltasMultF 直接支持）,
# 并与"两个 DiracDelta 关联函数相加"的方案逐位等价（关联函数对谱密度线性）。
# 运行: julia --startup-file=no --project=env2nd second_order_check_2mode.jl
# =============================================================================

using TEMPO, ImpurityModelBase, LinearAlgebra, Printf

BLAS.set_num_threads(4)

# ------------------------------ 模型参数 ------------------------------
const X = ComplexF64[0 1; 1 0]
const Y = ComplexF64[0 -im; im 0]
const Z = ComplexF64[1 0; 0 -1]
const HS = .37X + .19Y + .23Z
const SV = [-.7, .9]
const SMAT = ComplexF64[SV[1] 0; 0 SV[2]]
const RHO0 = (Matrix{ComplexF64}(I, 2, 2) + .3X + .4Y + .5Z) / 2

const MODES = [(ω = 1.3, g = 0.35), (ω = 0.7, g = 0.22)]   # 两模: 非对称、非公度
const BETA = 1.7
const TF = 1.6

nbar(ω) = 1 / (exp(BETA * ω) - 1)
# 解析关联函数 C(t) = Σ_j g_j²[(n̄_j+1)e^{-iω_j t} + n̄_j e^{+iω_j t}]
cfunc(t) = sum(m.g^2 * ((nbar(m.ω) + 1) * exp(-im * m.ω * t) + nbar(m.ω) * exp(im * m.ω * t)) for m in MODES)

# 解析 η 核（矩形/三角积分的逐模之和）
function eta_analytic(dt::Float64, lag::Int)
    z = 0.0 + 0.0im
    for m in MODES
        a = m.ω
        q0 = (1 - im * a * dt - exp(-im * a * dt)) / a^2          # 时序三角
        qr = 2 * (1 - cos(a * dt)) / a^2                          # 完整矩形
        z += m.g^2 * (lag == 0 ? ((nbar(a) + 1) * q0 + nbar(a) * q0') :
                                  qr * ((nbar(a) + 1) * exp(-im * a * lag * dt) + nbar(a) * exp(im * a * lag * dt)))
    end
    return z
end
if_analytic(sp, sm, dt) = begin
    z = 0.0 + 0.0im
    n = length(sp)
    for j in 1:n, k in 1:j
        e = eta_analytic(dt, j - k)
        z -= (sp[j] - sm[j]) * (e * sp[k] - e' * sm[k])
    end
    exp(z)
end

# ------------------------ 双模全空间 ED 基准 ------------------------
# 基底排序: kron(系统, kron(模1, 模2)); 偏迹用显式指标循环（避免 reshape 陷阱）
function bath_operators(nb::Int)
    # nb: 每模 Fock 截断（两模相同）; 基底 kron(模1, 模2)
    b1 = diagm(1 => sqrt.(1.0:(nb-1)))                             # 模1 湮灭
    b2 = diagm(1 => sqrt.(1.0:(nb-1)))                             # 模2 湮灭
    a1 = kron(b1, Matrix(I, nb, nb))
    a2 = kron(Matrix(I, nb, nb), b2)
    hb1, hb2 = MODES[1].ω, MODES[2].ω
    hb = kron(Matrix(Diagonal(hb1 * (0.0:(nb-1)))), Matrix(I, nb, nb)) +
              kron(Matrix(I, nb, nb), Matrix(Diagonal(hb2 * (0.0:(nb-1)))))
    vop = MODES[1].g * (a1 + a1') + MODES[2].g * (a2 + a2')
    p1 = exp.(-BETA * MODES[1].ω * (0.0:(nb-1))); p1 ./= sum(p1)
    p2 = exp.(-BETA * MODES[2].ω * (0.0:(nb-1))); p2 ./= sum(p2)
    rb = kron(Diagonal(p1), Diagonal(p2))
    return Matrix(hb), vop, rb, nb * nb
end

# 偏迹: ρ_{ab} = Σ_i R[(a-1)*nb+i, (b-1)*nb+i]，显式循环（已验证版本）
function ptrace_bath(R::AbstractMatrix, nb::Int)
    out = zeros(ComplexF64, 2, 2)
    for a in 1:2, b in 1:2, i in 1:nb
        out[a, b] += R[(a-1)*nb + i, (b-1)*nb + i]
    end
    return out
end

function ed_exact(nbt::Int)
    hb, vop, rb, nb = bath_operators(nbt)
    H = kron(HS, Matrix(I, nb, nb)) + kron(Matrix(I, 2, 2), hb) + kron(SMAT, vop)
    r0 = kron(RHO0, Matrix(rb))
    ev, V = eigen(Hermitian(H))
    U = V * Diagonal(exp.(-im * TF * ev)) * V'
    return ptrace_bath(U * r0 * U', nb)
end

function ed_split(nbt::Int, dt::Float64, order::Symbol)
    hb, vop, rb, nb = bath_operators(nbt)
    A = kron(exp(-im * dt * HS), Matrix(I, nb, nb))
    B = exp(-im * dt * (kron(Matrix(I, 2, 2), hb) + kron(SMAT, vop)))
    if order === :first
        step = A * B
    else
        h = exp(-im * dt / 2 * kron(HS, Matrix(I, nb, nb)))
        step = h * B * h
    end
    n = round(Int, TF / dt)
    U = step^n
    r0 = kron(RHO0, Matrix(rb))
    return ptrace_bath(U * r0 * U', nb)
end

# --------------------------- TEMPO（一阶内核） ---------------------------
function open_rdm(mps)
    L = length(mps)
    v = ones(ComplexF64, 1)
    for j in L:-1:3
        v = dropdims(sum(mps[j], dims=2), dims=2) * v
    end
    r = zeros(ComplexF64, 2, 2)
    for a in 1:2, b in 1:2
        r[a, b] = only(mps[1][:, a, :] * mps[2][:, b, :] * v) * scaling(mps)^L
    end
    r
end

# 双模浴: 原生 DiscreteSpectrum 入口（bosonic_Δt 通过 DeltasMultF 直接支持离散谱）
function corr2mode(lat)
    ws = [m.ω for m in MODES]
    fs = [m.g^2 for m in MODES]
    p = sortperm(ws)                       # DiscreteSpectrum 要求频率有序
    spec = DiscreteSpectrum(ws[p], fs[p])
    return correlationfunction(bosonicbath(spec, β=BETA), lat)
end
# 交叉验证: 与"两个 DiracDelta 关联函数相加"（关联函数对谱密度的线性性）一致
function corr2mode_summed(lat)
    return sum(correlationfunction(bosonicbath(DiracDelta(m.ω, α=m.g^2), β=BETA), lat) for m in MODES)
end

function build_if(n::Int; D=256, eps=1e-14)
    dt = TF / n
    lat = ADTLattice(N=n, δt=dt, contour=:real)
    corr = corr2mode(lat)
    trc = truncdimcutoff(D=D, ϵ=eps)
    mif = hybriddynamics(lat, corr, AdditiveHyb(SV), trunc=trc)
    return lat, mif, trc
end

function contract(mif, lat, trc, rho)
    k = sysdynamics(lat, ImpurityHamiltonian(HS), trunc=trc)
    initialstate!(k, lat, rho)
    m = mult(k, mif, trunc=trc)
    return open_rdm(m)
end

# 误差阶拟合: err ≈ C δt^p
fit_order(dts, errs) = begin
    x = log.(dts); y = log.(errs); n = length(x)
    (n * sum(x .* y) - sum(x) * sum(y)) / (n * sum(x .^ 2) - sum(x)^2)
end

# (A) IF(MPS) 与双模解析 η 的逐元素对比
function if_vs_analytic(n::Int=3)
    dt = TF / n
    lat = ADTLattice(N=n, δt=dt, contour=:real)
    corr = corr2mode(lat)
    trc = truncdimcutoff(D=256, ϵ=1e-14)
    mif = hybriddynamics(lat, corr, AdditiveHyb(SV), trunc=trc)
    function amplitude(mps, inds)
        v = ones(ComplexF64, 1)
        for j in eachindex(mps.data)
            v = vec(transpose(v) * mps[j][:, inds[j], :])
        end
        only(v) * scaling(mps)^length(mps)
    end
    errmax = 0.0
    for p in 0:(2^n-1), mm in 0:(2^n-1)
        inds = ones(Int, length(lat))
        for j in 1:n
            inds[lat[BosonicContourIndex(j, :+)]] = 1 + ((p >> (n-j)) & 1)
            inds[lat[BosonicContourIndex(j, :-)]] = 1 + ((mm >> (n-j)) & 1)
        end
        sp = [SV[1+((p >> (n-j)) & 1)] for j in 1:n]
        sm = [SV[1+((mm >> (n-j)) & 1)] for j in 1:n]
        errmax = max(errmax, abs(amplitude(mif, inds) - if_analytic(sp, sm, dt)))
    end
    return errmax
end

# ------------------------------- 主流程 -------------------------------
function main()
    t0 = time()
    println("=" ^ 78)
    println("二阶 QuAPI/TEMPO 离散化验证 —— 双模玻色浴: ω,g = ", [(m.ω, m.g) for m in MODES])
    println("=" ^ 78)

    nbt = 12   # 每模 Fock 截断
    exact = ed_exact(nbt)
    println("ED 浴截断检查  ‖ρ(12,12) − ρ(16,16)‖_F = ", @sprintf("%.3e", norm(ed_exact(16) - exact)))
    println("全空间 Lie  (δt=0.2) vs ED 误差 = ", @sprintf("%.3e", norm(ed_split(nbt, 0.2, :first) - exact)))
    println("全空间 Strang(δt=0.2) vs ED 误差 = ", @sprintf("%.3e", norm(ed_split(nbt, 0.2, :second) - exact)))

    # --- (A) IF vs 双模解析 η ---
    println("(A) N=3:  max |IF_MPS − IF_解析(双模)| = ", @sprintf("%.3e", if_vs_analytic(3)))
    # 离散谱入口 vs DiracDelta 叠加的等价性
    let lat = ADTLattice(N=4, δt=TF/4, contour=:real)
        c1, c2 = corr2mode(lat), corr2mode_summed(lat)
        eq = maximum(norm(Matrix(branch(c1, b1, b2)) - Matrix(branch(c2, b1, b2))) for b1 in (:+, :-), b2 in (:+, :-))
        println("(A') DiscreteSpectrum 入口 vs DiracDelta 叠加: max block diff = ", @sprintf("%.3e", eq))
    end
    println()

    # --- (B)(C) 误差阶 + 恒等式 ---
    Ns = [8, 16, 32, 64, 128]
    dts = TF ./ Ns
    e1s, e2s, d1s, d2s = Float64[], Float64[], Float64[], Float64[]
    io = open(joinpath(@__DIR__, "second_order_check_2mode_results.txt"), "w")
    hdr = @sprintf("%5s %9s %6s %14s %14s %14s %14s", "N", "dt", "bond",
                   "err_1st_vs_ED", "err_2nd_vs_ED", "|1st-Lie|", "|2nd-Strang|")
    println(hdr); println(io, hdr)
    for n in Ns
        dt = TF / n
        lat, mif, trc = build_if(n)
        r1 = contract(mif, lat, trc, RHO0)
        u = exp(-im * dt / 2 * HS)
        r2 = u' * contract(mif, lat, trc, u * RHO0 * u') * u
        rf = ed_split(nbt, dt, :first)
        rs = ed_split(nbt, dt, :second)
        e1, e2, d1, d2 = norm(r1-exact), norm(r2-exact), norm(r1-rf), norm(r2-rs)
        push!(e1s, e1); push!(e2s, e2); push!(d1s, d1); push!(d2s, d2)
        row = @sprintf("%5d %9.5f %6d %14.6e %14.6e %14.6e %14.6e", n, dt, bonddim(mif), e1, e2, d1, d2)
        println(row); flush(io); println(io, row)
    end
    println()
    p1a, p2a = fit_order(dts, e1s), fit_order(dts, e2s)
    p1b, p2b = fit_order(dts[1:4], e1s[1:4]), fit_order(dts[1:4], e2s[1:4])
    for s in (@sprintf("(B) 一阶误差阶: p = %.4f (N=8..128), %.4f (N=8..64)", p1a, p1b),
              @sprintf("    二阶误差阶: p = %.4f (N=8..128), %.4f (N=8..64)", p2a, p2b),
              @sprintf("(C) 恒等式: ‖ρ_1st−Lie‖ ≤ %.2e,  ‖ρ_2nd−Strang‖ ≤ %.2e (仅压缩残差)", maximum(d1s), maximum(d2s)))
        println(s); println(io, s)
    end
    close(io)
    println(); println("结果已写入 second_order_check_2mode_results.txt; 总耗时 ", round(time()-t0, digits=1), " s")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
