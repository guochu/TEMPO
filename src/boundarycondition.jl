"""
    boundarycondition(x::ADT, lattice::AbstractADTLattice; kwargs...)

Apply the lattice boundary condition on a copy of the ADT (calling [`boundarycondition!`](@ref)) and return the modified ADT.
"""
boundarycondition(x::ADT, lattice::AbstractADTLattice; kwargs...) = boundarycondition!(copy(x), lattice; kwargs...)


"""
    boundarycondition!(x::ADT, lattice::ImagADTLattice1Order)

Apply the boundary condition in place on an imaginary-time ADT lattice: connect the lattice points at the two ends of the imaginary-time contour with the identity matrix.
"""
function boundarycondition!(x::ADT, lattice::ImagADTLattice1Order)
	pos1, pos2 = index(lattice, 1), index(lattice, lattice.k)
	d = lattice.d
	return apply!(ADTTerm((pos1, pos2), _eye(d)), x)
	# canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	# return x
end



# the initial state of the impurity is maximally-mixed by default
"""
    boundarycondition!(x::ADT, lattice::RealADTLattice1Order; ρ₀::VecOrMat=ones(lattice.d))

Apply the boundary condition in place on a real-time ADT lattice: connect the ends of the `:+`/`:-` branches and apply the initial density matrix `ρ₀`.
By default `ρ₀` is the maximally mixed state (a vector of ones).
"""
function boundarycondition!(x::ADT, lattice::RealADTLattice1Order; ρ₀::VecOrMat=ones(lattice.d))
	d = lattice.d
	Is = _eye(d)
	pos1, pos2 = index(lattice, lattice.k, branch=:+), index(lattice, lattice.k, branch=:-)
	apply!(ADTTerm((pos1, pos2), Is), x)
	# canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return initialstate!(x, lattice, ρ₀)
end


"""
    initialstate!(x::ADT, lattice::RealADTLattice1Order, v0::AbstractVector)

Apply the initial state `v0` in vector form (diagonal placeholder) in place to the real-time ADT lattice; internally converts it to a diagonal density matrix and calls the matrix version.
"""
function initialstate!(x::ADT, lattice::RealADTLattice1Order, v0::AbstractVector)
	d = length(v0)
	m = zeros(eltype(v0), d, d)
	for i in 1:d
		m[i, i] = v0[i]
	end
	return initialstate!(x, lattice, m)
end
"""
    initialstate!(x::ADT, lattice::RealADTLattice1Order, ρ0::AbstractMatrix)

Apply the initial density matrix `ρ0` (normalized to `ρ0/tr(ρ0)`) in place to the first `:+`/`:-` index pair of the real-time ADT lattice.
"""
function initialstate!(x::ADT, lattice::RealADTLattice1Order, ρ0::AbstractMatrix)
	(size(ρ0, 1) == size(ρ0, 2) == lattice.d) || throw(DimensionMismatch("diagonal element size mismatch with phydim"))
	pos1, pos2 = index(lattice, 1, branch=:+), index(lattice, 1, branch=:-)
	return apply!(ADTTerm((pos1, pos2), ρ0/tr(ρ0)), x)
	# canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	# return x	
end


"""
    initialstate!(x::ProcessTensor, lattice::RealPTLattice1Order, ρ0::AbstractMatrix)

Apply the initial density matrix `ρ0` in place to the first `:+`/`:-` index pair of the real-time PT lattice (absorbing it into the process tensor via tensor decomposition).
"""
function initialstate!(x::ProcessTensor, lattice::RealPTLattice1Order, ρ0::AbstractMatrix)
	(size(ρ0, 1) == size(ρ0, 2) == lattice.d) || throw(DimensionMismatch("diagonal element size mismatch with phydim"))
	pos1, pos2 = index(lattice, 1, branch=:+), index(lattice, 1, branch=:-)
	@assert pos1 + 1 == pos2
	@tensor tmp[3,4,6,7] := ρ0[1,2] * x[pos1][3,4,5,1] * x[pos2][5,6,7,2] 
	u, s, v = tsvd!(tmp, (1,2), (3,4), trunc=DefaultIntegrationTruncation)

	I2 = one(ρ0)
	# I2 ./= sqrt(size(ρ0, 1))
	# I2 ./= tr(I2)
	s2 = Matrix(Diagonal(s))
	@tensor a[1,2,4,5,6] := u[1,2,3] * s2[3,4] * I2[5,6]
	@tensor b[1,4,2,3,5] := v[1,2,3] * I2[4,5]
	x[pos1] = tie(a, (1,1,2,1))
	x[pos2] = tie(b, (2,1,1,1))
	# canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	return x
end


"""
    boundarycondition!(x::ADT, lattice::MixedADTLattice1Order)

Apply the boundary condition in place on a mixed ADT lattice: connect the `:-`/`:+` ends and the two ends of the imaginary-time branch with the identity matrix.
"""
function boundarycondition!(x::ADT, lattice::MixedADTLattice1Order)
	d = lattice.d
	Is = _eye(d)

	pos1, pos2 = index(lattice, lattice.kt, branch=:-), index(lattice, lattice.kt, branch=:+)
	apply!(ADTTerm((pos1, pos2), Is), x)

	pos1, pos2 = index(lattice, 1, branch=:τ), index(lattice, 1, branch=:-)
	apply!(ADTTerm((pos1, pos2), Is), x)
	
	pos1, pos2 = index(lattice, 1, branch=:+), index(lattice, lattice.kτ, branch=:τ)
	return apply!(ADTTerm((pos1, pos2), Is), x)
	# canonicalize!(x, alg=Orthogonalize(trunc=trunc))
	# return x
end