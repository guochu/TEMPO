# =============================================================================
# 二阶 QuAPI/TEMPO 离散化验证 —— PT 路线 + 非对角耦合（单模玻色浴）
#
# 模型: H = H_S⊗1 + 1⊗ωa†a + S⊗g(a+a†), S 为工作基(σz 基)下的非对角厄米矩阵
#       S = [[0.2, 0.8], [0.8, -0.35]],  [H_S, S]≠0
#
# PT 工作流: PTLattice → hybriddynamics(NonAdditiveHyb(S)) → sysdynamics → mult
#            → initialstate! → rdm
#
# 已由 probe 确认: PT 路线 = Lie AB(左值) 排布, 非对角 IF 是分裂序列的精确表示
# (对全空间 ED-Lie 收敛到 ~4e-14)。二阶算法与 ADT 路线相同:
#     u = exp(-iδt H_S / 2);  ρ_in = u ρ0 u†;  运行一阶 PT;  ρ_out = u† ρ_out u
#
# 验证: (A) 一阶 ≈ 全空间 Lie, 二阶 ≈ 全空间 Strang (恒等式)
#       (B) 固定 T=1.6, 逐次减半 δt: 一阶 ~ δt¹, 二阶 ~ δt²
# 运行: julia --startup-file=no --project=env2nd second_order_check_nondiag.jl
# =============================================================================

using TEMPO, ImpurityModelBase, LinearAlgebra, Printf

BLAS.set_num_threads(4)

# ------------------------------ 模型参数 ------------------------------
const X = ComplexF64[0 1; 1 0]
const Y = ComplexF64[0 -im; im 0]
const Z = ComplexF64[1 0; 0 -1]
const HS = .37X + .19Y + .23Z
const SMAT = ComplexF64[0.2 0.8; 0.8 -0.35]      # 非对角厄米耦合算符
const RHO0 = (Matrix{ComplexF64}(I, 2, 2) + .3X + .4Y + .5Z) / 2
const OMEGA, G, BETA, TF = 1.3, 0.35, 1.7, 1.6

# ------------------------ 全空间 ED 基准（独立实现） ------------------------
function bath_matrices(nb::Int)
    b = diagm(1 => sqrt.(1.0:(nb-1)))
    hb = OMEGA * Matrix(Diagonal(0.0:(nb-1)))
    bb = G * (b + b')
    p = exp.(-BETA * OMEGA * (0.0:(nb-1))); p ./= sum(p)
    return hb, bb, Diagonal(p)
end
function ptrace_bath(R::AbstractMatrix, nb::Int)   # 显式循环, 已验证版本
    out = zeros(ComplexF64, 2, 2)
    for a in 1:2, c in 1:2, i in 1:nb
        out[a, c] += R[(a-1)*nb + i, (c-1)*nb + i]
    end
    out
end
function ed_exact(nb::Int)
    hb, bb, rb = bath_matrices(nb)
    H = kron(HS, Matrix(I, nb, nb)) + kron(Matrix(I, 2, 2), hb) + kron(SMAT, Matrix(bb))
    ev, V = eigen(Hermitian(H))
    U = V * Diagonal(exp.(-im * TF * ev)) * V'
    ptrace_bath(U * kron(RHO0, Matrix(rb)) * U', nb)
end
function ed_split(nb::Int, dt::Float64, order::Symbol)
    hb, bb, rb = bath_matrices(nb)
    A = kron(exp(-im * dt * HS), Matrix(I, nb, nb))
    B = exp(-im * dt * (kron(Matrix(I, 2, 2), hb) + kron(SMAT, Matrix(bb))))
    step = order === :first ? A * B :
           exp(-im * dt / 2 * kron(HS, Matrix(I, nb, nb))) * B * exp(-im * dt / 2 * kron(HS, Matrix(I, nb, nb)))
    U = step^round(Int, TF / dt)
    ptrace_bath(U * kron(RHO0, Matrix(rb)) * U', nb)
end

# --------------------------- PT 路线（一阶内核） ---------------------------
function run_pt(n::Int, rho; D=256, eps=1e-14)
    dt = TF / n
    lat = PTLattice(N=n, δt=dt, contour=:real)
    bath = bosonicbath(DiracDelta(OMEGA, α=G^2), β=BETA)
    corr = correlationfunction(bath, lat)
    trc = truncdimcutoff(D=D, ϵ=eps)
    mi = hybriddynamics(lat, corr, NonAdditiveHyb(SMAT), trunc=trc)
    k = sysdynamics(lat, ImpurityHamiltonian(HS), trunc=trc)
    m = mult(k, mi, trunc=trc)
    initialstate!(m, lat, rho)
    return rdm(lat, m), bonddim(mi)
end

# 误差阶拟合: err ≈ C δt^p
fit_order(dts, errs) = begin
    x = log.(dts); y = log.(errs); n = length(x)
    (n * sum(x .* y) - sum(x) * sum(y)) / (n * sum(x .^ 2) - sum(x)^2)
end

# ------------------------------- 主流程 -------------------------------
function main()
    t0 = time()
    println("=" ^ 78)
    println("二阶 QuAPI/TEMPO 离散化验证 —— PT 路线 + 非对角耦合 (NonAdditiveHyb)")
    println("S = [0.2 0.8; 0.8 -0.35] (非对角厄米), 单模浴 ω=$OMEGA, g=$G, β=$BETA, T=$TF")
    println("=" ^ 78)

    nb = 24
    exact = ed_exact(nb)
    println("ED 浴截断检查  ‖ρ(nb=24) − ρ(nb=32)‖_F = ", @sprintf("%.3e", norm(ed_exact(32) - exact)))
    println("全空间 Lie  (δt=0.2) vs ED 误差 = ", @sprintf("%.3e", norm(ed_split(nb, 0.2, :first) - exact)))
    println("全空间 Strang(δt=0.2) vs ED 误差 = ", @sprintf("%.3e", norm(ed_split(nb, 0.2, :second) - exact)))
    println()

    Ns = [8, 16, 32, 64, 128]
    dts = TF ./ Ns
    e1s, e2s, d1s, d2s = Float64[], Float64[], Float64[], Float64[]
    io = open(joinpath(@__DIR__, "second_order_check_nondiag_results.txt"), "w")
    hdr = @sprintf("%5s %9s %6s %14s %14s %14s %14s", "N", "dt", "bond",
                   "err_1st_vs_ED", "err_2nd_vs_ED", "|1st-Lie|", "|2nd-Strang|")
    println(hdr); println(io, hdr)
    for n in Ns
        dt = TF / n
        r1, bond = run_pt(n, RHO0)                     # 一阶: 包默认
        u = exp(-im * dt / 2 * HS)                     # 二阶: 端点半步旋转
        r2, _ = run_pt(n, u * RHO0 * u')
        r2 = u' * r2 * u
        rf = ed_split(nb, dt, :first)
        rs = ed_split(nb, dt, :second)
        e1, e2, d1, d2 = norm(r1-exact), norm(r2-exact), norm(r1-rf), norm(r2-rs)
        push!(e1s, e1); push!(e2s, e2); push!(d1s, d1); push!(d2s, d2)
        row = @sprintf("%5d %9.5f %6d %14.6e %14.6e %14.6e %14.6e", n, dt, bond, e1, e2, d1, d2)
        println(row); flush(io); println(io, row)
    end
    println()
    p1a, p2a = fit_order(dts, e1s), fit_order(dts, e2s)
    p1b, p2b = fit_order(dts[1:4], e1s[1:4]), fit_order(dts[1:4], e2s[1:4])
    for s in (@sprintf("(A) 恒等式: ‖ρ_1st−Lie‖ ≤ %.2e,  ‖ρ_2nd−Strang‖ ≤ %.2e (仅压缩残差)", maximum(d1s), maximum(d2s)),
              @sprintf("(B) 一阶误差阶: p = %.4f (N=8..128), %.4f (N=8..64)", p1a, p1b),
              @sprintf("    二阶误差阶: p = %.4f (N=8..128), %.4f (N=8..64)", p2a, p2b))
        println(s); println(io, s)
    end
    close(io)
    println(); println("结果已写入 second_order_check_nondiag_results.txt; 总耗时 ", round(time()-t0, digits=1), " s")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
