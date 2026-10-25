# 三合一检验: 对角耦合 z = diag(-0.7, 0.9) 走 NonAdditiveHyb
#   (a) partialif 路线 (应严格: 物理上=AdditiveHyb, ADT 已证 1e-13)
#   (b) XTRGIF 路线 k 扫描 (检查我之前 k=3 的偏差是否为微分步长未收敛)
#   (c) k 收敛后复检 σx 与 σ⁺ (手征 Wick 判定)
# 运行: julia --startup-file=no --project=env2nd probe_k_convergence.jl
using TEMPO, ImpurityModelBase, LinearAlgebra, Printf

const X = ComplexF64[0 1; 1 0]
const Y = ComplexF64[0 -im; im 0]
const Z = ComplexF64[1 0; 0 -1]
const HS = .37X + .19Y + .23Z
const ZD = ComplexF64[-.7 0; 0 .9]           # 对角耦合 (与 AdditiveHyb([-0.7,0.9]) 等价)
const SP = ComplexF64[0 1; 0 0]
const RHO0 = (Matrix{ComplexF64}(I, 2, 2) + .3X + .4Y + .5Z) / 2
const OMEGA, G, BETA = 1.3, 0.35, 1.7

ptrace(R, nb) = [sum(reshape(R, nb, 2, nb, 2)[i,a,i,b] for i in 1:nb) for a in 1:2, b in 1:2]
# op 厄米: H_int = G·op⊗(b+b†);  op=σ⁺: H_int = G(σ⁺⊗b + σ⁻⊗b†) (JC/NonDiagonalHyb)
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

contract(mi, lat, trc) = begin
    kk = sysdynamics(lat, ImpurityHamiltonian(HS), trunc=trc)
    m = mult(kk, mi, trunc=trc)
    initialstate!(m, lat, RHO0)
    rdm(lat, m)
end

function main()
    N, T = 4, 1.2
    lat = PTLattice(N=N, δt=T/N, d=2, contour=:real)
    corr = correlationfunction(bosonicbath(DiracDelta(OMEGA, α=G^2), β=BETA), lat)
    trc = truncdimcutoff(D=256, ϵ=1e-13)
    algexpan = OverDeterminedProny(n=20, tol=1e-10)
    D1 = DMRG1(trunc=trc, maxiter=10)

    rZ = ed_lie(ZD, N, T)
    println("== (a) partialif 路线, NonAdditiveHyb(z对角) ==")
    mi = hybriddynamics(lat, corr, NonAdditiveHyb(ZD), trunc=trc)
    @printf("   ‖PT − Lie‖ = %.3e   (bond=%d)\n", norm(contract(mi, lat, trc) - rZ), bonddim(mi)); flush(stdout)

    println("== (b) XTRGIF 路线 k 扫描, 同一 NonAdditiveHyb(z对角) ==")
    for k in (3, 5, 7)
        mi = hybriddynamics(lat, corr, NonAdditiveHyb(ZD),
                            XTRGIF(k=k, fast=true, alg_mult=D1, alg_expan=algexpan))
        @printf("   k=%d: ‖PT − Lie‖ = %.3e   (bond=%d)\n", k, norm(contract(mi, lat, trc) - rZ), bonddim(mi))
        flush(stdout)
    end

    println("== (c) k=7 复检 σx (NonAdditiveHyb) 与 σ⁺ (NonDiagonalHyb) ==")
    for (name, hyb) in (("σx", NonAdditiveHyb(X)), ("σ⁺", NonDiagonalHyb(SP)))
        r = ed_lie(name == "σx" ? X : SP, N, T)
        mi = hybriddynamics(lat, corr, hyb, XTRGIF(k=7, fast=true, alg_mult=D1, alg_expan=algexpan))
        @printf("   %s (k=7): ‖PT − Lie‖ = %.3e\n", name, norm(contract(mi, lat, trc) - r)); flush(stdout)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
