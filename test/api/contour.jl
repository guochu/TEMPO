# Contour indices

@testset "ContourIndex               " begin
	c = ContourIndex(3, branch=:+)
	@test c.j == 3
	@test branch(c) == :+
	@test ContourIndex(3, :+) == c
	@test ContourIndex(3, branch=:+) == c
	τ = ContourIndex(2)
	@test branch(τ) == :τ
	@test branch(ContourIndex(1, :-)) == :-
	# scalartype of arrays / numbers / lattice types
	@test scalartype(1.0) == Float64
	@test scalartype(randn(ComplexF64, 2, 2)) == ComplexF64
end
