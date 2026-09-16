# Fock orderings, layout/time-ordering styles and contour lattices (ADT and PT)

@testset "Fock orderings and styles" begin
	# ordering type hierarchy
	@test M2M1() isa FockOrdering
	@test M2M1() isa ImagFockOrdering
	@test M2m2M1m1() isa RealFockOrdering
	@test M2M1_m1M1m2M2() isa MixedFockOrdering
	# layout style
	@test LayoutStyle(M2M1()) isa TimeLocalLayout
	@test LayoutStyle(M2m2M1m1()) isa TimeLocalLayout
	@test LayoutStyle(M2M1_m1M1m2M2()) isa TimeLocalLayout
	# time-ordering styles
	@test ImaginaryTimeOrderingStyle(M2M1()) isa TimeDscending
	@test RealTimeOrderingStyle(M2m2M1m1()) isa TimeDscending
	@test TimeOrderingStyle(M2M1()) isa TimeDscending
	@test TimeOrderingStyle(M2m2M1m1()) isa TimeDscending
	@test RealTimeOrderingStyle(M2M1_m1M1m2M2()) isa TimeAscending
	@test ImaginaryTimeOrderingStyle(M2M1_m1M1m2M2()) isa TimeDscending
end

for (name, LatticeC) in (("ADT", ADTLattice), ("PT", PTLattice))
	local Lat = LatticeC
	@testset "$name lattice: imaginary time" begin
		lattice = Lat(N=2, δτ=0.1, d=3, contour=:imag, ordering=M2M1())
		@test LayoutStyle(lattice) isa TimeLocalLayout
		@test OrderingStyle(lattice) == M2M1()
		@test TimeOrderingStyle(lattice) isa TimeDscending
		@test lattice isa (name == "ADT" ? AbstractADTLattice : AbstractPTLattice)
		if name == "ADT"
			@test lattice isa ImagADTLattice
			@test lattice isa ImagADTLattice1Order
		else
			@test lattice isa ImagPTLattice
			@test lattice isa ImagPTLattice1Order
		end
		@test scalartype(lattice) == Float64
		@test length(lattice) == (name == "ADT" ? 3 : 2)
		@test lattice.d == 3
		@test lattice.N == 2
		@test lattice.β == 0.2
		@test lattice.δτ == 0.1
		@test lattice.τs == 0:0.1:0.2
		@test lattice.T == 5
		# index conventions: boundary site is the last lattice position
		for i in 1:length(lattice)
			@test index(lattice, i) == length(lattice) - i + 1
		end

		mps = vacuumstate(scalartype(lattice), lattice)
		@test scalartype(mps) == Float64
		if name == "ADT"
			@test integrate(mps) ≈ lattice.d^length(lattice) atol = 1.0e-6
		else
			@test integrate(lattice, mps) ≈ lattice.d atol = 1.0e-6
		end
	end

	@testset "$name lattice: real time" begin
		lattice = Lat(N=1, δt=0.1, contour=:real, ordering=M2m2M1m1())
		@test LayoutStyle(lattice) isa TimeLocalLayout
		@test OrderingStyle(lattice) == M2m2M1m1()
		@test TimeOrderingStyle(lattice) isa TimeDscending
		if name == "ADT"
			@test lattice isa RealADTLattice
			@test lattice isa RealADTLattice1Order
		else
			@test lattice isa RealPTLattice
			@test lattice isa RealPTLattice1Order
		end
		@test scalartype(lattice) == ComplexF64
		@test lattice.d == 2
		@test lattice.N == 1
		@test lattice.t == 0.1
		@test lattice.δt == 0.1
		@test lattice.ts == 0:0.1:0.1
		@test branches(lattice) == (:+, :-)
		# the ADT lattice includes the boundary step (k = N + 1), the PT lattice does not
		k = name == "ADT" ? lattice.N + 1 : lattice.N
		for i in 1:k
			if name == "ADT"
				# the boundary step (i = k) is the first lattice position
				@test index(lattice, i, branch=:+) == 2 * (k - i) + 1
				@test index(lattice, i, branch=:-) == 2 * (k - i) + 2
			else
				@test index(lattice, i, branch=:+) == 2 * i - 1
				@test index(lattice, i, branch=:-) == 2 * i
			end
		end

		mps = vacuumstate(scalartype(lattice), lattice)
		if name == "ADT"
			@test integrate(mps) ≈ lattice.d^length(lattice) atol = 1.0e-6
		else
			@test integrate(lattice, mps) ≈ lattice.d atol = 1.0e-6
		end
	end

	@testset "$name lattice: mixed time" begin
		lattice = Lat(Nt=1, δt=0.05, Nτ=2, δτ=0.1, contour=:mixed, ordering=M2M1_m1M1m2M2())
		@test LayoutStyle(lattice) isa TimeLocalLayout
		@test OrderingStyle(lattice) == M2M1_m1M1m2M2()
		# on the mixed contour the two branches have separate ordering styles
		@test ImaginaryTimeOrderingStyle(lattice) isa TimeDscending
		@test RealTimeOrderingStyle(lattice) isa TimeAscending
		if name == "ADT"
			@test lattice isa MixedADTLattice
			@test lattice isa MixedADTLattice1Order
		else
			@test lattice isa MixedPTLattice
			@test lattice isa MixedPTLattice1Order
		end
		@test scalartype(lattice) == ComplexF64
		@test lattice.Nt == 1
		@test lattice.Nτ == 2
		@test lattice.d == 2
		@test lattice.t == 0.05
		@test lattice.β == 0.2
		@test lattice.δt == 0.05
		@test lattice.ts == 0:0.05:0.05
		@test lattice.τs == 0:0.1:0.2
		@test branches(lattice) == (:+, :-, :τ)
		# imaginary time axis: the boundary point is the last τ position
		kτ = name == "ADT" ? lattice.Nτ + 1 : lattice.Nτ
		for i in 1:lattice.Nτ
			@test index(lattice, i, branch=:τ) == kτ + 1 - i
		end
		# real time axis
		@test index(lattice, 1, branch=:-) == kτ + 1
		@test index(lattice, 1, branch=:+) == kτ + 2

		mps = vacuumstate(scalartype(lattice), lattice)
		if name == "ADT"
			@test integrate(mps) ≈ lattice.d^length(lattice) atol = 1.0e-6
		else
			@test integrate(lattice, mps) ≈ lattice.d atol = 1.0e-6
		end
	end
end
