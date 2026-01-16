sysdynamics(gmps::ADT, lattice::RealADTLattice, model::ImpurityLindbladian, args...; kwargs...) = sysdynamics!(copy(gmps), lattice, model, args...; kwargs...)
function sysdynamics(lattice::RealADTLattice, model::ImpurityLindbladian, args...; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(model))
	sysdynamics!(vacuumstate(T, lattice), lattice, model, args...; kwargs...)
end 

function sysdynamics!(mps::ADT, lattice::RealADTLattice1Order, model::ImpurityLindbladian; trunc::TruncationScheme=DefaultKTruncation)
	alg = Orthogonalize(SVD(), trunc)
	U = _get_dissipative_adt_propagator(model.m, lattice.δt)
	for j in 1:lattice.N
		pos1, pos2 = index(lattice, j+1, branch=:+), index(lattice, j, branch=:+)
		pos3, pos4 = index(lattice, j, branch=:-), index(lattice, j+1, branch=:-)
		t = ADTTerm((pos1, pos2, pos3, pos4), U)
		apply!(t, mps)
		canonicalize!(mps, alg=alg)	
	end
	return mps
end

function sysdynamics!(mps::ADT, lattice::RealADTLattice1Order, model::ImpurityLindbladian, op::ContourOperator; trunc::TruncationScheme=DefaultKTruncation)
	alg = Orthogonalize(SVD(), trunc)
	U = _get_dissipative_adt_propagator(model.m, lattice.δt)
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

# single step
function sysdynamics!(mps::ADT, lattice::RealADTLattice1Order, model::ImpurityLindbladian, j::Int, op::AbstractMatrix; 
						trunc::TruncationScheme=DefaultKTruncation)
	alg = Orthogonalize(SVD(), trunc)
	U = _get_dissipative_adt_propagator(model.m, lattice.δt)
	pos1, pos2 = index(lattice, j+1, branch=:+), index(lattice, j, branch=:+)
	pos3, pos4 = index(lattice, j, branch=:-), index(lattice, j+1, branch=:-)
	U′ = _get_L_prime(U, j, op)
	t = ADTTerm((pos1, pos2, pos3, pos4), U′)
	apply!(t, mps)
	canonicalize!(mps, alg=alg)	
	return mps
end
function sysdynamics!(mps::ADT, lattice::RealADTLattice1Order, model::ImpurityLindbladian, j::Int; 
						trunc::TruncationScheme=DefaultKTruncation)
	alg = Orthogonalize(SVD(), trunc)
	U = _get_dissipative_adt_propagator(model.m, lattice.δt)
	pos1, pos2 = index(lattice, j+1, branch=:+), index(lattice, j, branch=:+)
	pos3, pos4 = index(lattice, j, branch=:-), index(lattice, j+1, branch=:-)
	t = ADTTerm((pos1, pos2, pos3, pos4), U)
	apply!(t, mps)
	canonicalize!(mps, alg=alg)	
	return mps
end


sysdynamics(gmps::ProcessTensor, lattice::RealPTLattice, model::ImpurityLindbladian; kwargs...) = sysdynamics!(copy(gmps), lattice, model; kwargs...)
function sysdynamics(lattice::RealPTLattice, model::ImpurityLindbladian; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(model))
	return sysdynamics!(vacuumstate(T, lattice), lattice, model; kwargs...)
end 


function sysdynamics!(mps::ProcessTensor, lattice::RealPTLattice1Order, model::ImpurityLindbladian; trunc::TruncationScheme=DefaultKTruncation)
	alg = Orthogonalize(SVD(), trunc)
	U = _get_dissipative_pt_propagator(model.m, lattice.δt)
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


# function _get_dissipative_adt_propagator(m4, δt)
# 	d = size(m4, 1)
# 	d2 = d * d
# 	m2 = reshape(m4, d2, d2)
# 	m2_exp = exp(m2 * δt)
# 	# return reshape(m2_exp, d,d,d,d)
# 	return permute(reshape(m2_exp, d,d,d,d), (1,3,4,2))
# end

function _get_dissipative_adt_propagator(m4, δt)
	r = _get_dissipative_pt_propagator(m4, δt)
	return permute(r, (1,3,4,2))
end

function _get_dissipative_pt_propagator(m4, δt)
	d = size(m4, 1)
	d2 = d * d
	m2 = reshape(m4, d2, d2)
	m2_exp = exp(m2 * δt)
	return reshape(m2_exp, d,d,d,d)
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