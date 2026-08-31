#=
    jctype.jl — JC-type spin-boson model (off-diagonal coupling, RWA)

Tutorial example: real-time dynamics of a spin-1/2 impurity coupled to a
bosonic bath through a *non-diagonal* (conjugate-pair) coupling, obtained from
the spin-boson model by the rotating-wave approximation (RWA). Demonstrates the
Process-Tensor (PT) TEMPO path of the package for off-diagonal system-bath
coupling on the real-time (Keldysh) contour.

Model (Jaynes-Cummings type):
    H = Δ σ_z/2 ⊗ 1 + 1 ⊗ H_bath + g Σ_k (σ₊ b_k + σ₋ b_k†)

with sub-Ohmic spectral density
    J(ω) = 2π α ω^s ω_c^{1-s},   0 ≤ ω ≤ ω_c,
inverse temperature β = 2.5, and initial impurity state ρ_imp = |+⟩⟨+|.

Pipeline:
  1. Build the real-time PT lattice: `PTLattice(N=Nt, δt=δt, contour=:real)`.
  2. Bath: `bath = bosonicbath(spectrum, β)`, `corr = correlationfunction(bath, lattice)`.
  3. Non-diagonal coupling: `hyb = NonDiagonalHyb(sp)` with sp = σ₊/2.
  4. Influence functional (translation-invariant, XTRG-style):
     `alg = TranslationInvariantIF(k=k, fast=true, algmult=DMRGMult1(...), algexpan=OverDeterminedProny(...))`,
     `mpsI = hybriddynamics(lattice, corr, hyb, alg)` (cached to `data/jc_... .mps`).
  5. System propagator: `mpsK = sysdynamics(lattice, model, trunc=trunc)` with
     `model = ImpurityHamiltonian(Δ .* z)`.
  6. Full tensor: `mps = mult!(mpsK, mpsI, trunc=trunc)`;
     `cache = environments(lattice, mps, ρ₀=ρimp)`.
  7. Observables: ⟨σ_x(t)⟩ at every time step via `ContourOperator` and
     `expectationvalue(m, cache)`.

Output:
  - `data/jc_... .mps`   : cached influence functional (ProcessTensor/MPO).
  - `result/jc_... .json`: `{"ts": [...], "z": [[Re,Im], ...]}`,
     real/imag parts of ⟨σ_x(t)⟩ on the time grid.

Run:   include("jctype.jl");  main(1.0)
=#

# push!(LOAD_PATH, "../../src")

using TEMPO, ImpurityModelBase, LinearAlgebra
using JSON, Serialization

# the spin-boson model after the rotating-wave approximation


function split_complexarray(a::Vector)
	L = length(a)
	r = similar(a, Float64, L, 2)
	r[:, 1] = real(a)
	r[:, 2] = imag(a)
	return r
end

function spin_half_matrices()
	s_SP = Array{Float64, 2}([0 0; 1 0])
	s_SM = Array{Float64, 2}([0 1; 0 0])
	s_Z = Array{Float64, 2}([-1 0; 0 1])
	s_x = s_SP+s_SM
	s_y = -im*(s_SP-s_SM)
	n = Array{Float64, 2}([0 0; 0 1])
	return Dict("x"=>s_x, "y"=>s_y, "z"=>s_Z, "+"=>s_SP, "-"=>s_SM, "n"=>n)
end

subomic_spectrum(w, α, s, wc) = 2π*α * w^s * wc^(1-s)


"""
    main(t; δt=0.05, Δ=1., β=2.5, α=0.1, s=0.5, wc=5., n=20, k=7, chi=100)

Run the JC-type (RWA spin-boson, off-diagonal coupling) example on the real-time contour.

The system is a spin-1/2 with Hamiltonian `Δ σ_z/2`, coupled to a sub-Ohmic
bosonic bath through the raising operator `σ₊` (`NonDiagonalHyb`). The
influence functional is built with the translation-invariant (XTRG-style) path
`TranslationInvariantIF`: the differential influence functional of width
`dt/2^k` is constructed by exponential expansion and WII time evolution, then
squared `k` times by the fast tree-bipartition scheme.  The observable
`⟨σ_x(t)⟩` is measured at every time step.

# Arguments
- `t::Real`: total real-time evolution time.
- `δt::Real=0.05`: real-time step size.
- `Δ::Real=1.0`: energy splitting (system Hamiltonian `Δ σ_z/2`).
- `β::Real=2.5`: inverse temperature of the bath.
- `α::Real=0.1`: dimensionless coupling strength (prefactor of the spectral density).
- `s::Real=0.5`: spectral density exponent (sub-Ohmic for `0 < s < 1`).
- `wc::Real=5.0`: spectral cutoff frequency.
- `n::Int=20`: number of Prony expansion terms (`OverDeterminedProny(n=n, tol=1e-8)`).
- `k::Int=7`: number of squaring steps of the fast tree-bipartition scheme
  (differential influence functional width `dt/2^k`).
- `chi::Int=100`: bond dimension cutoff `D` of the SVD truncation
  (`truncdimcutoff(D=chi, ϵ=1e-14)`).

# Returns
A `Dict` with
- `"ts"::Vector{Float64}`: time grid `0:δt:t`;
- `"z"::Matrix{Float64}`: `Nt × 2` matrix whose columns are the real and
  imaginary parts of `⟨σ_x(t)⟩` at each time step.

The same dictionary is also written to
`result/jc_β{β}_t{t}_dt{δt}_Δ{Δ}_α{α}_s{s}_wc{wc}_k{k}_n{n}_chi{chi}.json`,
and the influence functional is cached to `data/jc_... .mps` (computed only on
first run, loaded afterwards).
"""
function main(t; δt = 0.05, Δ = 1., β = 2.5, α=0.1, s=0.5, wc = 5., n=20, k=7, chi = 100)
	Nt = round(Int, t/δt)
	println("t=", t, " δt=", δt, " Δ=", Δ, " β=", β, " α=", α, " s=", s, " wc=", wc,  " n=", n, " k=", k, " chi=", chi)

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-14, add_back=0)
	# trunc2 = truncdimcutoff(D=2*chi, ϵ=1.0e-12, add_back=0)

	lattice = PTLattice(N=Nt, δt=δt, contour=:real)	
	println("numer of sites ", length(lattice))

	p = spin_half_matrices()
	x, y, z, sp = p["x"], p["y"], p["z"], p["+"]
	x, y, z, sp = 0.5 .* x, 0.5 .* y, 0.5 .* z, 0.5 .* sp

	mpspath = "data/jc_beta$(β)_t$(t)_dt$(δt)_alpha$(α)_s$(s)_wc$(wc)_k$(k)_n$(n)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		hyb = NonDiagonalHyb(sp)	
		spec = spectrum(w -> subomic_spectrum(w, α, s, wc), lb=0, ub=wc)
		bath = bosonicbath(spec, β=β)
		corr = correlationfunction(bath, lattice)
		algmult = DMRGMult1(trunc, maxiter=10)
		algexpan = OverDeterminedProny(n=n, tol=1.0e-8, verbosity=2)
		alg = TranslationInvariantIF(k=k, fast=true, algmult=algmult, algexpan=algexpan, verbosity=2)
		@time mpsI = hybriddynamics(lattice, corr, hyb, alg)
		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, mpsI)
	end

	println("mpsI bond dimension ", bond_dimension(mpsI))


	ρimp = [0.5 0.5;0.5 0.5]

	model = ImpurityHamiltonian(Δ .* z)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	@time mps = mult!(mpsK, mpsI, trunc=trunc)

	# tmp = initialstate!(copy(mps), lattice, ρimp)
	# Zval = integrate(lattice, tmp)
	cache = environments(lattice, mps, ρ₀ = ρimp)


	# local observables
	obs = ComplexF64[]
	for i in 1:Nt
		ind = ContourIndex(i, :+)
		m = ContourOperator(ind, x)
		# mps2 = apply!(m, lattice, copy(mps))
		# mps2 = initialstate!(mps2, lattice, ρimp)
		# v = integrate(lattice, mps2) / Zval
		v = expectationvalue(m, cache)
		push!(obs, v)
	end

	ts = [(i-1)*δt for i in 1:Nt]

	results = Dict("ts"=>ts, "z"=>split_complexarray(obs))

	data_path = "result/jc_beta$(β)_t$(t)_dt$(δt)_Delta$(Δ)_alpha$(α)_s$(s)_wc$(wc)_k$(k)_n$(n)_chi$(chi).json"

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return results	

end