push!(LOAD_PATH, "../../src")

using TEMPO, ImpurityModelBase, LinearAlgebra
using JSON, Serialization

# the standard spin-boson model

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

# const alphas = [0.04, 0.08, 0.12, 0.16]

function main(t; δt = 0.05, Δ = 1., β = 2.5, α=0.1, s=0.5, wc = 5., chi = 100)
	Nt = round(Int, t/δt)
	println("t=", t, " δt=", δt, " Δ=", Δ, " β=", β, " α=", α, " s=", s, " wc=", wc, " chi=", chi)

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-12, add_back=0)
	# trunc2 = truncdimcutoff(D=2*chi, ϵ=1.0e-12, add_back=0)

	lattice = ADTLattice(N=Nt, δt=δt, contour=:real)	
	println("numer of sites ", length(lattice))

	p = spin_half_matrices()
	x, y, z = p["x"], p["y"], p["z"]
	x, y, z = 0.5 .* x, 0.5 .* y, 0.5 .* z

	zdiag = [z[i,i] for i in 1:size(z, 1)]


	# read the mpsI if it has been already computed
	mpspath = "data/rabi_adt_beta$(β)_t$(t)_dt$(δt)_alpha$(α)_s$(s)_wc$(wc)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		hyb = AdditiveHyb(zdiag)	
		spec = spectrum(w -> subomic_spectrum(w, α, s, wc), lb=0, ub=wc)
		bath = bosonicbath(spec, β=β)
		corr = correlationfunction(bath, lattice)
		@time mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)
		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, mpsI)
	end

	println("mpsI bond dimension ", bond_dimension(mpsI))
	
	ρimp = [0 0; 0 1.]

	model = ImpurityHamiltonian(Δ .* x)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
	# @time mps = mult!(mpsK, mpsI, trunc=trunc2)
	# Zval = integrate(mps)
	cache = environments(lattice, mpsK, mpsI)


	# local observables, measure the local magnetization
	obs = ComplexF64[]
	for i in 1:Nt
		pos = index(lattice, i, branch=:+)
		m = ADTTerm(pos, zdiag)
		# mps2 = apply!(m, copy(mps))
		# v = integrate(mps2) / Zval
		v = expectationvalue(m, cache)
		push!(obs, v)
	end

	ts = [(i-1)*δt for i in 1:Nt]

	results = Dict("ts"=>ts, "z"=>split_complexarray(obs))

	data_path = "result/rabi_adt_beta$(β)_t$(t)_dt$(δt)_Delta$(Δ)_alpha$(α)_s$(s)_wc$(wc)_chi$(chi).json"

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return results	
end
