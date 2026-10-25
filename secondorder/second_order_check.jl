# =============================================================================
# 二阶 QuAPI/TEMPO 离散化验证（ADT 路线，对角耦合，时间无关 H_S）
#
# 二阶算法 = 端点半步旋转包装（IF 与所有内部系统整步原样复用）:
#     u     = exp(-i δt H_S / 2)
#     ρ_in  = u ρ0 u†        (输入端: 系统正向半步)
#     运行现有的一阶 TEMPO 演化
#     ρ_out = u† ρ_out u     (输出端: 系统反向半步)
#
# 验证判据:
#   (A) 固定 T=1.6, 逐次减半 δt:
#       一阶误差 ‖ρ-ρ_ED‖_F ~ δt^1,  二阶误差 ~ δt^2
#   (B) 一阶结果 ≈ 全空间 Lie 分裂, 二阶结果 ≈ 全空间 Strang 分裂
#       （恒等式 (uBu)^N = u†(AB)^N u 的直接检验，差值只来自 MPS 压缩）
#
# 模型（与 secondorder/ 比较报告一致）:
#   H = H_S ⊗ 1 + 1 ⊗ ω a†a + S ⊗ g(a + a†),  S = diag(-0.7, 0.9) 对角,
#   [H_S, S] ≠ 0, 单玻色模热浴。
# 运行: julia --startup-file=no --project=env2nd second_order_check.jl
# =============================================================================

using TEMPO, ImpurityModelBase, LinearAlgebra, Printf

BLAS.set_num_threads(4)

# ------------------------------ 模型参数 ------------------------------
const X = ComplexF64[0 1; 1 0]
const Y = ComplexF64[0 -im; im 0]
const Z = ComplexF64[1 0; 0 -1]
const HS = .37X + .19Y + .23Z              # 系统哈密顿量（与 S 不对易）
const SV = [-.7, .9]                       # 对角耦合算符本征值
const SMAT = ComplexF64[SV[1] 0; 0 SV[2]]
const RHO0 = (Matrix{ComplexF64}(I, 2, 2) + .3X + .4Y + .5Z) / 2
const OMEGA = 1.3                          # 模式频率
const G = 0.35                             # 耦合强度 (α = G²)
const BETA = 1.7
const TF = 1.6                             # 固定总时间

# ------------------------ 全空间 ED 基准（独立实现） ------------------------
function bath_matrices(nb::Int)
    b = diagm(1 => sqrt.(1.0:(nb-1)))      # 湮灭算符
    hb = OMEGA * Diagonal(0.0:(nb-1))      # 浴哈密顿量
    bb = G * (b + b')                      # 浴耦合算符 g(a + a†)
    p = exp.(-BETA * OMEGA * (0.0:(nb-1)))
    p ./= sum(p)                           # 热态
    return Matrix(hb), bb, Diagonal(p)
end

# 对浴偏迹: ρ_{ab} = Σ_i R[(a,i),(b,i)]。
# 注意 Julia 列主序: kron(系统, 浴) 的正确视图是 reshape(R, nb, 2, nb, 2)，
# 索引顺序为 (浴i, 系统a, 浴j, 系统b) —— 不能照抄 Python 的 reshape(2,nb,2,nb)（行主序）。
function ptrace_bath(R::AbstractMatrix, nb::Int)
    R4 = reshape(R, nb, 2, nb, 2)
    return [sum(R4[i, a, i, b] for i in 1:nb) for a in 1:2, b in 1:2]
end

# 精确传播（全空间对角化）
function ed_exact(nb::Int)
    hb, bb, rb = bath_matrices(nb)
    H = kron(HS, I(nb)) + kron(I(2), hb) + kron(SMAT, bb)
    r0 = kron(RHO0, Matrix(rb))
    ev, V = eigen(Hermitian(H))
    U = V * Diagonal(exp.(-im * TF * ev)) * V'
    return ptrace_bath(U * r0 * U', nb)
end

# 全空间 Lie / Strang 分裂（独立基准，用于检验恒等式）
function ed_split(nb::Int, dt::Float64, order::Symbol)
    hb, bb, rb = bath_matrices(nb)
    A = kron(exp(-im * dt * HS), I(nb))                    # 系统整步
    B = exp(-im * dt * (kron(I(2), hb) + kron(SMAT, bb)))  # 浴+耦合整步
    if order === :first
        step = A * B                                       # Lie
    else
        h = exp(-im * dt / 2 * kron(HS, I(nb)))
        step = h * B * h                                   # Strang
    end
    n = round(Int, TF / dt)
    U = step^n
    r0 = kron(RHO0, Matrix(rb))
    return ptrace_bath(U * r0 * U', nb)
end

# --------------------------- TEMPO（一阶内核） ---------------------------
# 末端 open pair 读出: + 支为行 (ket), − 支为列 (bra)
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
    return r
end

# 只构造一次 IF（每个 N 复用于一阶与二阶两种端点处理）
function build_if(n::Int; D=256, eps=1e-14)
    dt = TF / n
    lat = ADTLattice(N=n, δt=dt, contour=:real)
    bath = bosonicbath(DiracDelta(OMEGA, α=G^2), β=BETA)
    corr = correlationfunction(bath, lat)
    trc = truncdimcutoff(D=D, ϵ=eps)
    mif = hybriddynamics(lat, corr, AdditiveHyb(SV), trunc=trc)
    return lat, mif, trc
end

# 一阶收缩: 系统动力学 + 初态 + 与 IF 相乘 + open pair 读出
function contract(mif, lat, trc, rho)
    k = sysdynamics(lat, ImpurityHamiltonian(HS), trunc=trc)
    initialstate!(k, lat, rho)             # 末端保持 open（不加 trace δ）
    m = mult(k, mif, trunc=trc)
    return open_rdm(m)
end

# ====== 二阶算法：端点半步旋转（本脚本的核心实现） ======
function secondorder_wrap(rho_raw, u)
    return u' * rho_raw * u
end

# --------------------------- 误差阶拟合 ---------------------------
# err ≈ C δt^p  ⇒  log err = p log δt + c
function fit_order(dts, errs)
    x = log.(dts)
    y = log.(errs)
    n = length(x)
    return (n * sum(x .* y) - sum(x) * sum(y)) / (n * sum(x .^ 2) - sum(x)^2)
end

# ------------------------------- 主流程 -------------------------------
function main()
    t0 = time()
    println("=" ^ 78)
    println("二阶 QuAPI/TEMPO 离散化验证: 端点半步旋转 vs 一阶默认 (ADT 对角耦合)")
    println("=" ^ 78)

    # --- ED 基准与浴截断收敛 ---
    nb = 32
    exact = ed_exact(nb)
    println("模型: 2 能级 + 单玻色模; S = diag($SV); [H_S,S]≠0; T=$TF, β=$BETA, ω=$OMEGA, g=$G")
    println("ED 浴截断检查  ‖ρ(nb=32) − ρ(nb=48)‖_F = ", @sprintf("%.3e", norm(ed_exact(48) - exact)))
    println("全空间 Lie  (δt=0.2) vs ED 误差 = ", @sprintf("%.3e", norm(ed_split(nb, 0.2, :first) - exact)))
    println("全空间 Strang(δt=0.2) vs ED 误差 = ", @sprintf("%.3e", norm(ed_split(nb, 0.2, :second) - exact)))
    println()

    Ns = [8, 16, 32, 64, 128]
    dts = TF ./ Ns
    e1s = Float64[]      # 一阶 vs ED
    e2s = Float64[]      # 二阶 vs ED
    d1s = Float64[]      # 一阶 vs 全空间 Lie
    d2s = Float64[]      # 二阶 vs 全空间 Strang
    io = open(joinpath(@__DIR__, "second_order_check_results.txt"), "w")
    hdr = @sprintf("%5s %9s %6s %14s %14s %14s %14s %10s %10s",
                   "N", "dt", "bond", "err_1st_vs_ED", "err_2nd_vs_ED",
                   "|1st-Lie|", "|2nd-Strang|", "|tr1-1|", "|tr2-1|")
    println(hdr)
    println(io, hdr)

    for n in Ns
        dt = TF / n
        lat, mif, trc = build_if(n)
        bond = bonddim(mif)

        # 一阶：包默认
        r1 = contract(mif, lat, trc, RHO0)

        # 二阶：同一 IF，端点半步旋转
        u = exp(-im * dt / 2 * HS)
        r2 = secondorder_wrap(contract(mif, lat, trc, u * RHO0 * u'), u)

        # 独立基准
        rf = ed_split(nb, dt, :first)
        rs = ed_split(nb, dt, :second)

        e1 = norm(r1 - exact); e2 = norm(r2 - exact)
        d1 = norm(r1 - rf);    d2 = norm(r2 - rs)
        push!(e1s, e1); push!(e2s, e2); push!(d1s, d1); push!(d2s, d2)
        row = @sprintf("%5d %9.5f %6d %14.6e %14.6e %14.6e %14.6e %10.1e %10.1e",
                       n, dt, bond, e1, e2, d1, d2, abs(tr(r1) - 1), abs(tr(r2) - 1))
        println(row); flush(io)
        println(io, row)
    end

    println()
    # 阶数拟合（全部点，以及去掉最细网格 N=128——压缩残差可能在那里开始主导）
    p1a = fit_order(dts, e1s);       p2a = fit_order(dts, e2s)
    p1b = fit_order(dts[1:4], e1s[1:4]); p2b = fit_order(dts[1:4], e2s[1:4])
    msg1 = @sprintf("一阶误差阶:  p = %.4f (N=8..128),  %.4f (N=8..64)", p1a, p1b)
    msg2 = @sprintf("二阶误差阶:  p = %.4f (N=8..128),  %.4f (N=8..64)", p2a, p2b)
    msg3 = @sprintf("恒等式检验:  ‖ρ_1st − Lie‖_F   ≤ %.2e  (仅 MPS 压缩残差)", maximum(d1s))
    msg4 = @sprintf("            ‖ρ_2nd − Strang‖_F ≤ %.2e  (仅 MPS 压缩残差)", maximum(d2s))
    for s in (msg1, msg2, msg3, msg4)
        println(s); println(io, s)
    end
    close(io)
    println()
    println("结果已写入 second_order_check_results.txt; 总耗时 ", round(time() - t0, digits=1), " s")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
