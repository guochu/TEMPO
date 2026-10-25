# 最终判定: 固定 k=7, 扫 N — σ⁺(NonDiagonalHyb) 的偏差斜率
#   斜率 2 ⇒ 纯微分截断, ti_localop 等时修正充分 (TEMPO NonDiagonalHyb 严格)
#   斜率 1 ⇒ O(δt) Wick 残留 (手征耦合不严格)
# 对照组: σx (NonAdditiveHyb, 已知严格) 应为斜率 2
# 运行: julia --startup-file=no --project=env2nd probe_wick_final.jl
using TEMPO, ImpurityModelBase, LinearAlgebra, Printf

const X = ComplexF64[0 1; 1 0]
const Y = ComplexF64[0 -im; im 0]
const Z = ComplexF64[1 0; 0 -1]
const HS = .37X + .19Y + .23Z
const SP = ComplexF64[0 1; 0 0]
const RHO0 = (Matrix{ComplexF64}(I, 2, 2) + .3X + .4Y + .5Z) / 2
const OMEGA, G, BETA = 1.3, 0.35, 1.7

ptrace(R, nb) = [sum(reshape(R, nb, 2, nb, 2)[i,a,i,b] for i in 1:nb) for a in 1:2, b in 1:2]
function ed_lie(op, N, T; nb=20)
    h = T / N
    bth = diagm(1 => sqrt.(1.0:(nb-1)))
    hb = OMEGA * Matrix(Diagonal(0.0:(nb-1)))
    V = ishermitian(op) ? G * kron(op, bth + bth') :
                         G * (kron(op, bth) + kron(op', bth'))
    rb = exp.(-BETA * OMEGA * (0.0:(nb-1))); rb ./= sum(rb)
    A = kron(exp(-im*h*HS), Matrix(I, nb, nb))
    B = exp(-im*h*(kron(Matrix(I,2,2), hb) + V))
    U = (A * B)^N
    ptrace(U * kron(RHO0, Diagonal(rb)) * U', nb)
end

function dev(hyb, op, N, T, trc, D1, algexpan)
    lat = PTLattice(N=N, δt=T/N, d=2, contour=:real)
    corr = correlationfunction(bosonicbath(DiracDelta(OMEGA, α=G^2), β=BETA), lat)
    mi = hybriddynamics(lat, corr, hyb, XTRGIF(k=7, fast=true, alg_mult=D1, alg_expan=algexpan))
    kk = sysdynamics(lat, ImpurityHamiltonian(HS), trunc=trc)
    m = mult(kk, mi, trunc=trc)
    initialstate!(m, lat, RHO0)
    norm(rdm(lat, m) - ed_lie(op, N, T))
end

function main()
    T = 1.2
    trc = truncdimcutoff(D=256, ϵ=1e-13)
    D1 = DMRG1(trunc=trc, maxiter=10)
    algexpan = OverDeterminedProny(n=20, tol=1e-10)
    for (name, hyb, op) in (("σx", NonAdditiveHyb(X), X), ("σ⁺", NonDiagonalHyb(SP), SP))
        errs = Float64[]
        for N in (4, 8, 16)
            e = dev(hyb, op, N, T, trc, D1, algexpan)
            push!(errs, e)
            @printf("%s  N=%2d  δt=%.4f:  ‖PT − Lie‖ = %.3e\n", name, N, T/N, e); flush(stdout)
        end
        p = (log(errs[3]) - log(errs[1])) / (log(T/16) - log(T/4))
        @printf("%s  拟合斜率 (N=4→16): p = %.2f   [2 ⇒ 无 Wick 残留; 1 ⇒ O(δt) Wick]\n\n", name, p); flush(stdout)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
