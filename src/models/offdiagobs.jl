


function _sysdynamics_util!(gmps::ADT, lattice::AbstractADTLattice, model::BosonicImpurity, branch::Symbol, 
							N::Int, cts::ContourOperator; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	U = propagator(model, lattice, branch)
	# data = decompose_to_mps(U)
	alg = Orthogonalize(SVD(), trunc)
	for j in 1:N
		a, b = (branch == :-) ? (j, j+1) : (j+1, j)
        pos1, pos2 = index(lattice, a, branch=branch), index(lattice, b, branch=branch)
        U′ = _get_U_prime(lattice, U, j, branch, cts)
        t = ADTTerm((pos1, pos2), U′)
        apply!(t, gmps)
        canonicalize!(gmps, alg=alg)			
	end
	return gmps
end

function _get_U_prime(lattice, U, j::Int, bh::Symbol, cts::ContourOperator)
	idx, ops = cts.indices, cts.ops
	a = ContourIndex(j, bh)
	k = findfirst(x->x==a, idx)

	U′ = U
	if !isnothing(k)
		if bh == :-
			U′ = ops[k] * U
		else
			U′ = U * ops[k]  
		end
	end

	return U′
end




