println("------------------------------------")
println("|         Quantum Control          |")
println("------------------------------------")


function random_unitary(d)
	x = randn(ComplexF64, d, d)
	x = x + x'
	return exp(im .* x)
end


@testset "Single spin" begin
	Ω = 0.5
	N = 10
	δt = 0.5
	t = N * δt
	β = 1
	chi = 50
	tol = 1.0e-6
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)


	ρimp = _rand_dm(2)
	seq = [random_unitary(2) for i in 1:N]
	
	# hop = Ω .* [0 1; 1 0]
	hop = Ω .* pauli_y()
	model = ImpurityHamiltonian(hop)
	Uop = exp((-im .* δt) .* hop)

	ρ = ρimp
	for i in 1:N
		lattice = PTLattice(N = i, δt=δt, contour=:real)
		adt = sysdynamics(lattice, model, trunc=trunc)
		inds = [[ContourIndex(j, :+) for j in 1:i]; [ContourIndex(j, :-) for j in 1:i]]
		ddseq = seq[1:i]
		# ddseq2 = [Matrix(item') for item in ddseq]
		op = ContourOperator(inds, [ddseq; adjoint.(ddseq)]) 
		apply!(op, lattice, adt)
		initialstate!(adt, lattice, ρimp)
		rho1 = rdm(lattice, adt)
		rho1 ./= 2


		rho2 = seq[i] * ρ * seq[i]'
		rho2 = Uop * rho2 * Uop'

		err = distance(rho1, rho2) / norm(rho2)
		# println("i=", i, " error=", err)
		@test err < tol
		ρ = rho2
	end

end

@testset "Rabi model" begin

	Ω = 0.5
	N = 15
	δt = 0.05
	β = 1
	t = N * δt
	chi = 100
	d = 50
	tol = 1.0e-2
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)
	

	p = spin_half_matrices()
	x, y, z = p["x"], p["y"], p["z"]
	hop = Ω .* z
	Is = one(x)
	Ib = one(zeros(d, d))
	model = ImpurityHamiltonian(hop)

	Hbarebath = bosondensityoperator(d=d)
	a = bosonaoperator(d=d)
	H = kron(hop, Ib) + kron(Is, Hbarebath) + kron(y, a' + a)

	Uop = exp((-im .* δt) .* H)
	ρbath = exp(-β .* Hbarebath)
	ρimp = _rand_dm(2)
	ρ = kron(ρimp, ρbath) 
	ρ ./= tr(ρ)

	seq = [random_unitary(2) for i in 1:N]
	hyb = NonAdditiveHyb(y)

	spec = DiracDelta(1)
	bath = bosonicbath(spec, β=β)

	for i in 1:N
		lattice = PTLattice(N = i, δt=δt, contour=:real)
		mpsK = sysdynamics(lattice, model, trunc=trunc)
		corr = correlationfunction(bath, lattice)
		mpsI = hybriddynamics_naive(lattice, corr, hyb, trunc=trunc)
		adt = mult!(mpsK, mpsI, trunc=trunc)
		inds = [[ContourIndex(j, :+) for j in 1:i]; [ContourIndex(j, :-) for j in 1:i]]
		ddseq = seq[1:i]
		op = ContourOperator(inds, [ddseq; adjoint.(ddseq)]) 
		apply!(op, lattice, adt)
		initialstate!(adt, lattice, ρimp)
		rho1 = rdm(lattice, adt)
		rho1 ./= 2

		ρ4 = reshape(ρ, d, 2, d, 2)
		@tensor tmp4[3,1,4,6] := seq[i][1,2] * ρ4[3,2,4,5] * conj(seq[i][6,5])
		# tmp2 = seq[i] * ρ * seq[i]'
		tmp2 = reshape(tmp4, d*2, d*2)
		ρ = Uop * tmp2 * Uop'
		tmp4 = reshape(ρ, d, 2, d, 2)
		@tensor rho2[2,3] := tmp4[1,2,1,3]

		# println(rho1)
		# println(rho2)
		err = distance(rho1, rho2) / norm(rho2)
		# println("i=", i, " error=", err)
		@test err < tol
		
	end
end
