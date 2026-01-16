# ADT
"""
	sysdynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::ImpurityHamiltonian; kwargs...)

The inplace version of the function sysdynamics, which applys the bare impurity onto the given GMPS
"""
function sysdynamics!(gmps::ADT, lattice::ImagADTLattice, model::ImpurityHamiltonian, args...; trunc::TruncationScheme=DefaultKTruncation)
	return sysdynamics_imaginary!(gmps, lattice, model, args...; trunc=trunc)
end 


function sysdynamics!(gmps::ADT, lattice::RealADTLattice, model::ImpurityHamiltonian, args...; 
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		sysdynamics_forward!(gmps, lattice, model, args...; trunc=trunc)
		return sysdynamics_backward!(gmps, lattice, model, args...; trunc=trunc)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		return (branch == :+) ? sysdynamics_forward!(gmps, lattice, model, args...; trunc=trunc) : sysdynamics_backward!(gmps, lattice, model, args...; trunc=trunc)
	end
end 

function sysdynamics!(gmps::ADT, lattice::MixedADTLattice, model::ImpurityHamiltonian, args...; 
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		sysdynamics_forward!(gmps, lattice, model, args...; trunc=trunc)
		sysdynamics_backward!(gmps, lattice, model, args...; trunc=trunc)
		return sysdynamics_imaginary!(gmps, lattice, model, args...; trunc=trunc)
	else
		if branch == :+
			return sysdynamics_forward!(gmps, lattice, model, args...; trunc=trunc)
		elseif branch == :-
			return sysdynamics_backward!(gmps, lattice, model, args...; trunc=trunc)
		else
			(branch == :τ) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
			return sysdynamics_imaginary!(gmps, lattice, model, args...; trunc=trunc)
		end
	end
end 

sysdynamics_forward!(mps::ADT, lattice::AbstractADTLattice, model::ImpurityHamiltonian, args...; trunc::TruncationScheme=DefaultKTruncation) = _sysdynamics_util!(
						mps, lattice, model, :+, lattice.Nt, args...; trunc=trunc)
sysdynamics_backward!(mps::ADT, lattice::AbstractADTLattice, model::ImpurityHamiltonian, args...; trunc::TruncationScheme=DefaultKTruncation) = _sysdynamics_util!(
						mps, lattice, model, :-, lattice.Nt, args...; trunc=trunc)
sysdynamics_imaginary!(mps::ADT, lattice::AbstractADTLattice, model::ImpurityHamiltonian, args...; trunc::TruncationScheme=DefaultKTruncation) = _sysdynamics_util!(
						mps, lattice, model, :τ, lattice.Nτ, args...; trunc=trunc)

function sysdynamics!(mps::ADT, lattice::AbstractADTLattice, model::ImpurityHamiltonian, ind::ContourIndex, op::AbstractMatrix; 
						trunc::TruncationScheme=DefaultKTruncation)
	bh = branch(ind)
	U = propagator(model, lattice, bh)
	alg = Orthogonalize(SVD(), trunc)
	if bh == :-
		a, b = j, j+1
		U′ = op * U
	else
		a, b = j+1, j
		U′ = U * op
	end
	pos1, pos2 = index(lattice, a, branch=bh), index(lattice, b, branch=bh)
	t = ADTTerm((pos1, pos2), U′)
	apply!(t, mps)
	canonicalize!(mps, alg=alg)
	return mps
end
function sysdynamics!(mps::ADT, lattice::AbstractADTLattice, model::ImpurityHamiltonian, ind::ContourIndex; 
						trunc::TruncationScheme=DefaultKTruncation)
	bh = branch(ind)
	U = propagator(model, lattice, bh)
	alg = Orthogonalize(SVD(), trunc)
	if bh == :-
		a, b = j, j+1
	else
		a, b = j+1, j
	end
	pos1, pos2 = index(lattice, a, branch=bh), index(lattice, b, branch=bh)
	t = ADTTerm((pos1, pos2), U)
	apply!(t, mps)
	canonicalize!(mps, alg=alg)
	return mps
end


function _sysdynamics_util!(gmps::ADT, lattice::AbstractADTLattice, model::ImpurityHamiltonian, branch::Symbol, N::Int; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	U = propagator(model, lattice, branch)
	# data = decompose_to_mps(U)
	alg = Orthogonalize(SVD(), trunc)
	for j in 1:N
		a, b = (branch == :-) ? (j, j+1) : (j+1, j)
        pos1, pos2 = index(lattice, a, branch=branch), index(lattice, b, branch=branch)
        t = ADTTerm((pos1, pos2), U)
        apply!(t, gmps)
        canonicalize!(gmps, alg=alg)			
	end
	return gmps
end

function _sysdynamics_util!(gmps::ADT, lattice::AbstractADTLattice, model::ImpurityHamiltonian, branch::Symbol, 
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


function _get_propagator(h, lattice, b::Symbol)
	if b == :τ
		return exp(-lattice.δτ .* h)
	elseif b == :+
		return exp(-im*lattice.δt .* h)
	else
		return exp(im*lattice.δt .* h)
	end
end