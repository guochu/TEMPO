
# function gf(lattice::ImagADTLattice, model::AbstractImpurityOperator, op1::AbstractMatrix, op2::AbstractMatrix, mpsI::ADT; 
# 			b1::Symbol=:τ, b2::Symbol=:τ, trunc::TruncationScheme=DefaultKTruncation)
# 	T = promote_type(scalartype(lattice), scalartype(mpsI))
# 	mpsK = vaccumstate(T, lattice)
# 	rs = Vector{T}(undef, lattice.N)
# 	left = ones(T, 1)
# 	left′ = left
# 	for i in 1:lattice.N
# 		mpsK = sys
# 	end
# end