println("------------------------------------")
println("|       PT Zipup Integrate         |")
println("------------------------------------")


@testset "PT-zipup integrate: imaginary time" begin
	N = 4
	tol = 1.0e-6
	for d in (2,3)
		lattice = PTLattice(N=N, δτ=0.1, d=d, contour=:imag)
		L = length(lattice)
		for T in (Float64, ComplexF64)
			pt1 = randompt(T, L, d=d, D=3)
			canonicalize!(pt1)
			pt2 = randompt(T, L, d=d, D=4)
			canonicalize!(pt2)
			pt = pt1 * pt2
			v1 = integrate(lattice, pt)
			v2 = integrate(lattice, pt1, pt2)
			@test abs(v1 - v2) / abs(v1) < tol
		end
	end
end

@testset "PT-zipup integrate: real time" begin
	N = 4
	tol = 1.0e-6
	for d in (2,3)
		lattice = PTLattice(N=N, δt=0.1, d=d, contour=:real)
		L = length(lattice)
		pt1 = randompt(ComplexF64, L, d=d, D=3)
		canonicalize!(pt1)
		pt2 = randompt(ComplexF64, L, d=d, D=4)
		canonicalize!(pt2)
		pt = pt1 * pt2
		v1 = integrate(lattice, pt)
		v2 = integrate(lattice, pt1, pt2)
		@test abs(v1 - v2) / abs(v1) < tol
	end
end


@testset "PT-zipup integrate: mixed time" begin
	Nτ = 4
	Nt = 2
	tol = 1.0e-6
	for d in (2,3)
		lattice = PTLattice(Nt=Nt, Nτ=Nτ, δt=0.1, δτ=0.1, d=d, contour=:mixed)
		L = length(lattice)
		pt1 = randompt(ComplexF64, L, d=d, D=3)
		canonicalize!(pt1)
		pt2 = randompt(ComplexF64, L, d=d, D=4)
		canonicalize!(pt2)
		pt = pt1 * pt2
		v1 = integrate(lattice, pt)
		v2 = integrate(lattice, pt1, pt2)
		@test abs(v1 - v2) / abs(v1) < tol
	end
end