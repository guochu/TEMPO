# Contour indices

@testset "BosonicContourIndex" begin
	c = BosonicContourIndex(3, branch=:+)
	@test c.j == 3
	@test branch(c) == :+
	@test BosonicContourIndex(3, :+) == c
	@test BosonicContourIndex(3, branch=:+) == c
	τ = BosonicContourIndex(2)
	@test branch(τ) == :τ
	@test branch(BosonicContourIndex(1, :-)) == :-
	# scalartype of arrays / numbers / lattice types
	@test scalartype(1.0) == Float64
	@test scalartype(randn(ComplexF64, 2, 2)) == ComplexF64
end
