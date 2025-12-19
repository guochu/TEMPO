boundarycondition(x::ADT, lattice::AbstractADTLattice; kwargs...) = boundarycondition!(copy(x), lattice; kwargs...)


function boundarycondition!(x::ADT, lattice::ImagADTLattice1Order; trunc::TruncationScheme=DefaultIntegrationTruncation)
	pos1, pos2 = index(lattice, 1), index(lattice, lattice.k)
	d = lattice.d
	apply!(ADTTerm((pos1, pos2), _eye(d)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x
end



# the initial state of the impurity is maximally-mixed by default
function boundarycondition!(x::ADT, lattice::RealADTLattice1Order; ρ₀::VecOrMat=ones(lattice.d), trunc::TruncationScheme=DefaultIntegrationTruncation)
	d = lattice.d
	Is = _eye(d)
	pos1, pos2 = index(lattice, lattice.k, branch=:+), index(lattice, lattice.k, branch=:-)
	apply!(ADTTerm((pos1, pos2), Is), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return initialstate!(x, lattice, ρ₀, trunc=trunc)
end


function initialstate!(x::ADT, lattice::RealADTLattice1Order, v0::AbstractVector; kwargs...)
	d = length(v0)
	m = zeros(eltype(v0), d, d)
	for i in 1:d
		m[i, i] = v0[i]
	end
	return initialstate!(x, lattice, m; kwargs...)
end
function initialstate!(x::ADT, lattice::RealADTLattice1Order, ρ0::AbstractMatrix; trunc::TruncationScheme=DefaultIntegrationTruncation)
	(size(ρ0, 1) == size(ρ0, 2) == lattice.d) || throw(DimensionMismatch("diagonal element size mismatch with phydim"))
	pos1, pos2 = index(lattice, 1, branch=:+), index(lattice, 1, branch=:-)
	apply!(ADTTerm((pos1, pos2), ρ0/tr(ρ0)), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x	
end



function boundarycondition!(x::ADT, lattice::MixedADTLattice1Order; trunc::TruncationScheme=DefaultIntegrationTruncation)
	d = lattice.d
	Is = _eye(d)

	pos1, pos2 = index(lattice, lattice.kt, branch=:-), index(lattice, lattice.kt, branch=:+)
	apply!(ADTTerm((pos1, pos2), Is), x)

	pos1, pos2 = index(lattice, 1, branch=:τ), index(lattice, 1, branch=:-)
	apply!(ADTTerm((pos1, pos2), Is), x)
	
	pos1, pos2 = index(lattice, 1, branch=:+), index(lattice, lattice.kτ, branch=:τ)
	apply!(ADTTerm((pos1, pos2), Is), x)
	canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x
end