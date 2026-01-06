# dissipative impurity
struct DissipativeImpurity
	m::Array{ComplexF64, 4}
end

phydim(h::DissipativeImpurity) = size(h.m, 1)

DissipativeImpurity(d::Int) = DissipativeImpurity(zeros(ComplexF64, d, d, d, d))
TO.scalartype(::Type{DissipativeImpurity}) = ComplexF64
DissipativeImpurity(L::LindbladOperator) = DissipativeImpurity(L.m)
DissipativeImpurity(H::AbstractMatrix, jumpops::Vector{<:AbstractMatrix}) = DissipativeImpurity(lindbladoperator(H, jumpops))

# function lindbladoperator(H::AbstractMatrix, jumpops::Vector{<:AbstractMatrix})
# 	I2 = one(H)
# 	@tensor L[1,4,2,3] := -im * H[1,2] * I2[3,4]
# 	@tensor L[1,4,2,3] += im * I2[1,2] * H[3,4]
# 	for A in jumpops
# 		AdagA = A' * A
# 		@tensor L[1,4,2,3] -= AdagA[1,2] * I2[3,4]
# 		@tensor L[1,4,2,3] -= I2[1,2] * AdagA[3,4]
# 		@tensor L[1,4,2,3] += 2.0 * A[1,2] * A'[3,4]
# 	end
# 	return L
# end

# function check_lindbladoperator(H::AbstractMatrix, jumpops::Vector{<:AbstractMatrix}, rho::AbstractMatrix)
# 	d = size(H, 1)
# 	# rho = randn(ComplexF64, d, d)

# 	rho1 = -im * H * rho + im * rho * H
# 	for A in jumpops
# 		AdagA = A' * A
# 		rho1 -= (AdagA * rho + rho * AdagA)
# 		rho1 += 2 * A * rho * A'
# 	end
	
# 	L = lindbladoperator(H, jumpops)
# 	d2 = d * d
# 	rho2 = reshape(reshape(L, d2, d2) * reshape(rho, d2), d, d)
# 	return rho1, rho2
# end


sysdynamics(gmps::ADT, lattice::AbstractADTLattice, model::DissipativeImpurity, args...; kwargs...) = sysdynamics!(copy(gmps), lattice, model, args...; kwargs...)
function sysdynamics(lattice::AbstractADTLattice, model::DissipativeImpurity, args...; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(model))
	sysdynamics!(vacuumstate(T, lattice), lattice, model, args...; kwargs...)
end 

function sysdynamics!(mps::ADT, lattice::RealADTLattice1Order, model::DissipativeImpurity; trunc::TruncationScheme=DefaultKTruncation)
	alg = Orthogonalize(SVD(), trunc)
	U = _get_dissipative_propagator(model.m, lattice.δt)
	for j in 1:lattice.N
		pos1, pos2 = index(lattice, j+1, branch=:+), index(lattice, j, branch=:+)
		pos3, pos4 = index(lattice, j, branch=:-), index(lattice, j+1, branch=:-)
		t = ADTTerm((pos1, pos2, pos3, pos4), U)
		apply!(t, mps)
		canonicalize!(mps, alg=alg)	
	end
	return mps
end

function sysdynamics!(mps::ADT, lattice::RealADTLattice1Order, model::DissipativeImpurity, op::ContourOperator; trunc::TruncationScheme=DefaultKTruncation)
	alg = Orthogonalize(SVD(), trunc)
	U = _get_dissipative_propagator(model.m, lattice.δt)
	for j in 1:lattice.N
		pos1, pos2 = index(lattice, j+1, branch=:+), index(lattice, j, branch=:+)
		pos3, pos4 = index(lattice, j, branch=:-), index(lattice, j+1, branch=:-)
		U′ = _get_L_prime(U, j, op)
		t = ADTTerm((pos1, pos2, pos3, pos4), U′)
		apply!(t, mps)
		canonicalize!(mps, alg=alg)	
	end
	return mps
end


sysdynamics(gmps::ProcessTensor, lattice::AbstractPTLattice, model::DissipativeImpurity; kwargs...) = sysdynamics!(copy(gmps), lattice, model; kwargs...)
function sysdynamics(lattice::AbstractPTLattice, model::DissipativeImpurity; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(model))
	return sysdynamics!(vacuumstate(T, lattice), lattice, model; kwargs...)
end 


function sysdynamics!(mps::ProcessTensor, lattice::RealPTLattice1Order, model::DissipativeImpurity; trunc::TruncationScheme=DefaultKTruncation)
	alg = Orthogonalize(SVD(), trunc)
	U = _get_dissipative_propagator(model.m, lattice.δt)
	for j in 1:lattice.N
		ind1, ind2 = ContourIndex(j, :+), ContourIndex(j, :-)
		# t = ContourOperator((ind1, ind2), U)
		# apply!(t, lattice, mps)	
		t = FockTermS((lattice[ind1], lattice[ind2]), U)
		apply!(t, mps)	
	end
	canonicalize!(mps, alg=alg)
	return mps
end

function _get_dissipative_propagator(m4, δt)
	d = size(m4, 1)
	d2 = d * d
	m2 = reshape(m4, d2, d2)
	m2_exp = exp(m2 * δt)
	# return reshape(m2_exp, d,d,d,d)
	return permute(reshape(m2_exp, d,d,d,d), (1,3,4,2))
end

function _get_L_prime(U, j::Int, cts::ContourOperator)
	idx, ops = cts.indices, cts.ops
	U′ = U
	for bh in (:+, :-)
		a = ContourIndex(j, bh)
		k = findfirst(x->x==a, idx)

		if !isnothing(k)
			if bh == :-
				@tensor tmp[3,4,1,5] := ops[k][1,2] * U′[3,4,2,5]
			else
				@tensor tmp[1,5,3,4] := U′[1,2,3,4] * ops[k][2,5]  
			end
			U′ = tmp
		end
	end
	return U′
end