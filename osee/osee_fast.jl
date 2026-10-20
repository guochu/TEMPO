using TEMPO, ImpurityModelBase, FiniteMPSAlgorithms
using FiniteMPSAlgorithms: todense
using LinearAlgebra

omic_spectrum(w, α, wc) = 2α * w * exp(-w / wc)

N, α, wc = 4, 0.1, 2.0
δτ = 0.1; β = N * δτ
spec = spectrum(w -> omic_spectrum(w, α, wc), lb=0, ub=wc)
bath = bosonicbath(spec, β=β)

# Δ 与 Fock 截断 d 无关
lattice3 = PTLattice(N=N, δτ=δτ, d=3, contour=:imag)
corr = correlationfunction(bath, lattice3)
Δ = Matrix(corr.data)
println("Δ ="); display(round.(Δ, digits=8))
println("\neigvals(Δ) = ", round.(eigvals(Δ), digits=10))

# ---------- 精确 quadratic 模型：用 ImpurityModelBase 构造 H = Σ Δ_ij a†_i a_j ----------
function exact_osee(Δ, d)
    N = size(Δ, 1)
    terms = [adaga(i, j, coeff=Δ[i, j]) for i in 1:N for j in 1:N if Δ[i, j] != 0]
    h = quadratichamiltonian(N, terms)
    H = bosonoperator(h; d=d)
    F = exp(H)
    return [op_entropy(F, d, N, l) for l in 1:N-1], norm(F)
end

function op_entropy(F, d, N, l)
    R = reshape(F, ntuple(_ -> d, 2N)...)
    rowdims = [N - s + 1 for s in 1:N]
    coldims = [N + (N - s + 1) for s in 1:N]
    perm = vcat(vcat(rowdims[1:l], coldims[1:l]), vcat(rowdims[l+1:N], coldims[l+1:N]))
    R2 = permutedims(R, perm)
    Mmat = reshape(R2, d^(2l), d^(2(N - l)))
    sv = svdvals(Mmat); p = sv .^ 2; p ./= sum(p); p = p[p .> 1e-18]
    return -sum(p .* log.(p))
end

println("\n--- 精确 quadratic 模型 exp(+H)，d 扫描（d^N 维稠密） ---")
for d in (3, 4, 5, 6, 8)
    t = @elapsed (ss, nf) = exact_osee(Δ, d)
    println("d=$d  S = ", round.(ss, digits=10), "   ||F||=$(round(nf,digits=6))   ($(round(t,digits=2))s)")
end

# ---------- 快速解析公式：Δ 的 off-diagonal block SVD 主导阶 ----------
println("\n--- 快速公式：S_l ≈ 2 Σ_k s_k^2 (1 - ln s_k^2)，s_k = svdvals(Δ[L,R]) ---")
for l in 1:N-1
    blk = Δ[1:l, l+1:N]
    s = svdvals(blk)
    S = 2 * sum(x^2 * (1 - log(x^2)) for x in s if x > 0)
    println("l=$l  sv(block)=", round.(s, digits=8), "   S_fast = ", round(S, digits=10))
end

# ---------- TEMPO PT 在不同 Fock 截断 d 下 ----------
println("\n--- TEMPO PT (XTRGIF k=5, D=200) ---")
trunc = truncdimcutoff(D=200, ϵ=1e-14, add_back=0)
algmult = TEMPO.DMRG1(trunc=trunc)
alg_expan = OverDeterminedProny(n=3, tol=1e-12)
for d in (3, 6)
    lattice = PTLattice(N=N, δτ=δτ, d=d, contour=:imag)
    corr_d = correlationfunction(bath, lattice)
    alg = XTRGIF(k=5, fast=true, alg_mult=algmult, alg_expan=alg_expan)
    t = @elapsed pt = hybriddynamics(lattice, corr_d, NonDiagonalHyb(bosonadagoperator(d=d)), alg)
    Fpt = todense(pt.parent; order=:msb)
    ss = [op_entropy(Fpt, d, N, l) for l in 1:N-1]
    println("d=$d  bonddim=$(bonddim(pt))  PT bond S = ", round.(ss, digits=10), "   ($(round(t,digits=2))s)")
end
