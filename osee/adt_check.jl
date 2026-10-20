using TEMPO, ImpurityModelBase, LinearAlgebra

omic_spectrum(w, α, wc) = 2α * w * exp(-w / wc)
N, α, wc = 4, 0.1, 2.0
δτ = 0.1; β = N * δτ
d = 3

spec = spectrum(w -> omic_spectrum(w, α, wc), lb=0, ub=wc)
bath = bosonicbath(spec, β=β)

latpt = PTLattice(N=N, δτ=δτ, d=d, contour=:imag)
Δ = Matrix(correlationfunction(bath, latpt).data)
println("Δ =\n", round.(Δ, digits=8))

# ---------- 精确 dense：F = exp(sgn * Σ Δ_ij n_i n_j) ----------
nloc = Matrix(Diagonal(collect(0.0:d-1)))
function onsite(op, k)
	mats = [Matrix{Float64}(I, d, d) for _ in 1:N]
	mats[k] = op
	return kron(mats...)
end
Hd = zeros(d^N, d^N)
for i in 1:N, j in 1:N
	Δ[i, j] == 0 && continue
	Hd .+= Δ[i, j] .* (onsite(nloc, i) * onsite(nloc, j))
end

function op_entropy(F, l)
	R = reshape(F, ntuple(_ -> d, 2N)...)
	rowdims = [N - s + 1 for s in 1:N]
	coldims = [N + (N - s + 1) for s in 1:N]
	perm = vcat(vcat(rowdims[1:l], coldims[1:l]), vcat(rowdims[l+1:N], coldims[l+1:N]))
	R2 = permutedims(R, perm)
	Mmat = reshape(R2, d^(2l), d^(2(N - l)))
	sv = svdvals(Mmat); p = sv .^ 2; p ./= sum(p); p = p[p .> 1e-18]
	return -sum(p .* log.(p))
end

for sgn in (+1, -1)
	F = exp(sgn * Hd)
	println("exact exp($(sgn > 0 ? "+" : "-")H_diag)  S = ",
		round.([op_entropy(F, l) for l in 1:N-1], digits=10), "  norm=", round(norm(F), digits=6))
end

# ---------- TEMPO ADT ----------
lat = ADTLattice(N=N, δτ=δτ, d=d, contour=:imag)
corr = correlationfunction(bath, lat)
hyb = AdditiveHyb(collect(0.0:d-1))          # 耦合算符 = n̂
trunc = truncdimcutoff(D=200, ϵ=1e-14, add_back=0)
algmult = TEMPO.DMRG1(trunc=trunc)
alg_expan = OverDeterminedProny(n=3, tol=1e-12)
alg = XTRGIF(k=5, fast=true, alg_mult=algmult, alg_expan=alg_expan)

t = @elapsed adt = hybriddynamics(lat, corr, hyb, alg)
println("ADT build time = ", round(t, digits=2), " s, bonddim = ", bonddim(adt))
println("ADT physical dims = ", [size(x, 2) for x in adt.data])
Sadt = []
for (i, s) in enumerate(adt.s)
	sv = s ./ norm(s)
	p = sv .^ 2
	p = p[p .> 1e-18]
	push!(Sadt, -sum(p .* log.(p)))
end
println("ADT bond entropies = ", round.(Sadt, digits=10))
