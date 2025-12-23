println("------------------------------------")
println("|           PT PartialIF           |")
println("------------------------------------")

function _rand_hermitian(::Type{T}, d::Int) where {T<:Number}
	m = randn(T, d, d)
	return m + m'
end

@testset "PartialIF: imaginary time" begin
	δτ=0.1
	N = 2
	β = N * δτ
	tol = 1.0e-4

	spec = Leggett(d=1, ωc=1)
	bath = bosonicbath(spec, β=β)
	corr = Δτ(bath, N=N, δτ=δτ)

	for d in [2,]
		lattice = PTLattice(N=N, δτ=δτ, d=d, contour=:imag)	
		# op = _rand_hermitian(scalartype(lattice), d)
		op = [1 0; 0 -1]
		hyb = NonAdditiveHyb(op)
		for i in 1:lattice.N
			println("i = ", i)
			ind1 = ContourIndex(i)
			p1 = partialif(lattice, ind1, corr, hyb)
			p2 = partialif_naive(lattice, ind1, corr, hyb)

			println(distance(p1, p2), " ", norm(p1), " ", norm(p2))

			@test distance(p1, p2) / norm(p1) < tol
		end
	end

end


# @testset "PartialIF: real time" begin
# 	δt=0.1
# 	N = 10
# 	t = N * δt
# 	β = 2
# 	tol = 1.0e-4

# 	spec = Leggett(d=1, ωc=1)
# 	bath = bosonicbath(spec, β=β)
# 	corr = Δt(bath, N=N, t=t)

# 	for d in [2,3]
# 		lattice = ADTLattice(N=N, δt=δt, d=d, contour=:real)	
# 		op = rand(d)
# 		hyb = AdditiveHyb(op)
# 		for i in 1:lattice.N, b1 in branches(lattice)
# 			ind1 = ContourIndex(i, branch=b1)
# 			p1 = partialif(lattice, ind1, corr, hyb)
# 			p2 = partialif_naive(lattice, ind1, corr, hyb)

# 			@test distance(p1, p2) / norm(p1) < tol
# 		end
# 	end
# end


# @testset "PartialIF: mixed time" begin
# 	δτ=0.1
# 	Nτ = 10
# 	δt=0.05
# 	Nt = 5
# 	t = Nt * δt
# 	β = 2
# 	tol = 1.0e-4

# 	spec = Leggett(d=1, ωc=1)
# 	bath = bosonicbath(spec, β=β)
# 	corr = Δm(bath, Nτ=Nτ, Nt=Nt, t=t)

# 	for d in [2,3]
# 		lattice = ADTLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, d=d, contour=:mixed)	
# 		op = rand(d)
# 		hyb = AdditiveHyb(op)
# 		for b1 in branches(lattice)
# 			N = (b1 == :τ) ? lattice.Nτ : lattice.Nt
# 			for i in 1:N
# 				ind1 = ContourIndex(i, branch=b1)
# 				p1 = partialif(lattice, ind1, corr, hyb)
# 				p2 = partialif_naive(lattice, ind1, corr, hyb)

# 				@test distance(p1, p2) / norm(p1) < tol
# 			end
# 		end
# 	end
# end
