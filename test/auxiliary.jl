println("------------------------------------")
println("|            Auxiliary             |")
println("------------------------------------")

function prodmpo(::Type{T}, ds::Vector{Int}, positions::Vector{Int}, ops::Vector{<:AbstractMatrix}) where {T <: Number}
	(length(positions) == length(ops)) || throw(DimensionMismatch("positions and ops size mismatch"))
	(length(Set(positions)) == length(positions)) || throw(ArgumentError("multiple n̂ on the same position not allowed"))
	L = length(ds)
	mpotensors = Vector{Array{T, 4}}(undef, L)
	for i in 1:L
		pos = findfirst(x->x==i, positions)
		dj = ds[i]
		if isnothing(pos)
			mj = one(zeros(T, dj, dj))
		else
			mj = ops[pos]
		end
		mpotensors[i] = reshape(mj, (1,dj,1,dj))
	end
	return ProcessTensor(mpotensors)
end
prodmpo(ds::Vector{Int}, positions::Vector{Int}, ops::Vector{<:AbstractMatrix}) = prodmpo(Float64, ds, positions, ops)
function prodmpo(L::Int, positions::Vector{Int}, ops::Vector{<:AbstractMatrix})
	d = size(ops[1], 1)
	return prodmpo([d for i in 1:L], positions, ops)
end

function longrange_xxz(J, Jzz, hz, α, p)
	sp, sm, z = p["+"], p["-"], p["z"]
	C = [sp, sm, z]
	B = [2*J * sp', 2*J * sm', Jzz * z]
	terms = []
	for (a1, a2) in zip(C, B)
		push!(terms, ExponentialDecayTerm(a1, a2, α=exp(-α)))
	end
	return SchurMPOTensor(hz * z, [terms...])
end

function longrange_xxz_ham(L, hz, J, Jzz, α, p)
	sp, sm, z = p["+"], p["-"], p["z"]
	mpo = prodmpo(L, [1], [hz * z])
	for i in 2:L
		mpo += prodmpo(L, [i], [hz * z])
	end
	canonicalize!(mpo)
	for i in 1:L
	    for j in i+1:L
	    	coeff = exp(-α*(j-i))
	    	mpo += prodmpo(L, [i, j], [2*J*coeff*sp, sp'])
	    	mpo += prodmpo(L, [i, j], [2*J*coeff*sm, sm'])
	    	mpo += prodmpo(L, [i, j], [Jzz*coeff*z, z])
	    	canonicalize!(mpo)
	    end
	end
	return mpo
end

function longrange_xxz_mpoham(L, hz, J, Jzz, α, p)
	# the last term of J and Jzz not used
	mpo = MPOHamiltonian([longrange_xxz(J, Jzz, hz, α, p) for i in 1:L])
end

function powlaw_xxz(L, J, Jzz, hz, α, p)
	sp, sm, z = p["+"], p["-"], p["z"]
	C = [sp, sm]
	B = [2*J * sp', 2*J * sm']
	terms = []
	for (a1, a2) in zip(C, B)
		push!(terms, ExponentialDecayTerm(a1, a2, α=exp(-1)))
	end
	append!(terms, exponential_expansion(PowerlawDecayTerm(z, Jzz*z, α=α), len=L, alg=PronyExpansion(tol=1.0e-8)))
	return SchurMPOTensor(hz * z, [terms...])
end

powerlaw_xxz_mpoham(L, J, Jzz, hz, α, p) = MPOHamiltonian([powlaw_xxz(L, J, Jzz, hz, α, p) for i in 1:L])


function powerlaw_xxz_ham(L, J, Jzz, hz, α, p)
	sp, sm, z = p["+"], p["-"], p["z"]
	mpo = prodmpo(L, [1], [hz * z])
	for i in 2:L
		mpo += prodmpo(L, [i], [hz * z])
	end
	canonicalize!(mpo)
	for i in 1:L
	    for j in i+1:L
	    	coeff = exp(-(j-i))
	    	mpo += prodmpo(L, [i, j], [2*J*coeff*sp, sp'])
	    	mpo += prodmpo(L, [i, j], [2*J*coeff*sm, sm'])
	    	coeff = (j-i)^α
	    	mpo += prodmpo(L, [i, j], [Jzz*coeff*z, z])
	    	canonicalize!(mpo)
	    end
	end
	return mpo
end

@testset "Exponential expansion    " begin
	L = 100
	atol = 1.0e-5
	for alpha in (-2, -2.5, -3)
		xdata = [convert(Float64, i) for i in 1:L]
		ydata = [1.3 * x^alpha for x in xdata]
		xs1, lambdas1 = exponential_expansion(ydata, PronyExpansion(n=20,tol=atol))
		@test expansion_error(ydata, xs1, lambdas1) / norm(ydata) < atol
		xs2, lambdas2 = exponential_expansion(ydata, DeterminedPronyExpansion(n=20,tol=atol))
		@test expansion_error(ydata, xs2, lambdas2) / norm(ydata) < atol
		# xs3, lambdas3 = exponential_expansion(ydata, PronyExpansion2(n=20,atol=atol))
		# @test expansion_error(ydata, xs3, lambdas3) < atol
		# xs2, lambdas2 = exponential_expansion(ydata, LsqExpansion(atol=atol))
		# @test expansion_error(ydata, xs2, lambdas2) < atol
	end
	L = 500
	xdata = [convert(Float64, i) for i in 1:L]
	ydata = [1.3 * 0.7^x + 0.7 * 0.5^x - 1.1 * 0.8^x + 0.1*0.95^x for x in xdata]
	xs1, lambdas1 = exponential_expansion(ydata, PronyExpansion(n=20,tol=atol))
	@test expansion_error(ydata, xs1, lambdas1) / norm(ydata) < atol
	for stepsize in 2:6
		xs2, lambdas2 = exponential_expansion(ydata[1:stepsize:L], PronyExpansion(n=20, stepsize=stepsize, tol=atol))
		@test expansion_error(ydata, xs2, lambdas2) / norm(ydata) < atol
	end
	L = 500
	xdata = [convert(Float64, i) for i in 1:L]
	ydata = [(1.3+0.2im) * (0.7+0.3im)^x + (0.7+1.1im) * (0.5+0.1im)^x - (1.1-0.3im) * (0.8-0.5im)^x + (0.1-0.2)*(0.95+0.2im)^x for x in xdata]
	for alg in (PronyExpansion(n=20,tol=atol),)
		xs, lambdas = exponential_expansion(ydata, alg)
		@test expansion_error(ydata, xs, lambdas) / norm(ydata) < atol
	end
end

@testset "MPOHamiltonian: long-range XXZ        " begin
	p = spin_half_matrices()
	for L in (2, 3, 4)
		hz = 0.8
		J = 1
		Jzz = 1.2
		α = 0.9
		h1 = ProcessTensor(tompotensors(longrange_xxz_mpoham(L, hz, J, Jzz, α, p)))
		@test length(h1) == L
		h2 = longrange_xxz_ham(L, hz, J, Jzz, α, p)
		@test length(h2) == L
		@test distance(h1, h2) / norm(h2) < 1.0e-6	
	end
end

@testset "MPOHamiltonian: power-law XXZ      " begin
	p = spin_half_matrices()
	L = 20
	α = -2.5
	hz = 0.8
	J = 1
	Jzz = 1.2
	h1 = ProcessTensor(tompotensors(powerlaw_xxz_mpoham(L, hz, J, Jzz, α, p)))
	h2 = powerlaw_xxz_ham(L, hz, J, Jzz, α, p)

	@test distance(h1, h2) / norm(h2) < 1.0e-5
end
