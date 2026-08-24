

struct ADTExpectationCache{M<:ADT, G<:Tuple, L<:AbstractADTLattice, Hl, Hr} 
	A::M
	Bs::G
	lattice::L
	hleft::Hl
	hright::Hr	
end

Base.length(x::ADTExpectationCache) = length(x.lattice)
"""
	Zvalue(cache)

Return the partition function Z (the scalar obtained from contracting the left and right environments). For an `ADTExpectationCache` this is the rightmost left environment; for a `PTExpectationCache` it is the contraction of the left and right environments at the boundary.

See also [`expectationvalue`](@ref) and [`environments`](@ref).
"""
Zvalue(x::ADTExpectationCache) = only(x.hleft[end])
leftenv(x::ADTExpectationCache, j::Int) = x.hleft[j]
rightenv(x::ADTExpectationCache, j::Int) = x.hright[j+1]

function ADTExpectationCache(lattice::AbstractADTLattice, As::Tuple)
	(all(v->length(v)==length(lattice), As)) || throw(DimensionMismatch())

	L = length(lattice)

	left = l_LL(As...)
	hleft = Vector{typeof(left)}(undef, L+1)
	hleft[1] = left
	for i in 1:L
		left = left * TransferMatrix(i, As...)
		hleft[i+1] = left
	end
	right = r_RR(As...)
	hright = Vector{typeof(right)}(undef, L+1)
	hright[L+1] = right
	for i in L:-1:1
		right = TransferMatrix(i, As...) * right
		hright[i] = right
	end
	return ADTExpectationCache(first(As), Base.tail(As), lattice, hleft, hright)
end

"""
	environments(lattice, A, B...)
	environments(lattice::RealPTLattice, A::ProcessTensor, B::ProcessTensor...; ρ₀=_eye(phydim(lattice)))

Compute the left and right (boundary) environments required for expectation values, returning an expectation value cache.

ADT usage: `environments(lattice::AbstractADTLattice, A::ADT, B::ADT...)`, where `A` is the state being observed (if A and B coincide, this gives the partition function), returning an `ADTExpectationCache`.

PT usage: `environments(lattice::ImagPTLattice/MixedPTLattice, A::ProcessTensor, B::ProcessTensor...)` or `environments(lattice::RealPTLattice, A, B...; ρ₀=...)`. For real-time PT lattices, the keyword argument `ρ₀::AbstractMatrix` specifies the initial density matrix (default `_eye(phydim(lattice))`), used as the right boundary condition.

# Arguments
- `lattice`: contour lattice.
- `A`: the process tensor / augmented density tensor being observed.
- `B...`: additional tensors entering the contraction (usually influence functionals).
- `ρ₀::AbstractMatrix`: (real-time PT only) initial density matrix.

# Returns
An `ADTExpectationCache` or `PTExpectationCache`, for use with [`expectationvalue`](@ref), [`expectation`](@ref) and [`Zvalue`](@ref).
"""
environments(lattice::AbstractADTLattice, A::ADT, B::Vararg{ADT}) = ADTExpectationCache(lattice, (A, B...))

"""
	expectationvalue(m, cache)

Compute the expectation value of the operator `m`: ⟨m⟩ = expectation(m, cache) / Zvalue(cache). `m` can be an `ADTTerm` (ADT version) or a `ContourOperator`/`AbstractFockTerm` (PT version).

# Returns
The scalar expectation value.
"""
expectationvalue(m::ADTTerm, cache::ADTExpectationCache) = expectation(m, cache) / Zvalue(cache)
"""
	expectation(m, cache)

Compute the unnormalized expectation value (the numerator ⟨ψ̄|m|ψ⟩), usually used together with `Zvalue(cache)`. The ADT version accepts an `ADTTerm`; the PT version accepts an `AbstractFockTerm`.

See also [`expectationvalue`](@ref).
"""
function expectation(m::ADTTerm, cache::ADTExpectationCache)
	j, k = m.positions[1], m.positions[end]
	left = leftenv(cache, j)  
	right = rightenv(cache, k) 
	A2 = apply!(m, copy(cache.A))
	for tj in k:-1:j
		right = TransferMatrix(tj, A2, cache.Bs...) * right 
	end	
	# @tensor r = left[1] * right[1]
	# return r
	return contract_center(left, right)
end

# function contract_center(left::Vector, right::Vector)
# 	@tensor r = left[1] * right[1]
# 	return r
# end

# function contract_center(left::Matrix, right::Matrix)
# 	@tensor r = left[1, 2] * right[1, 2]
# 	return r
# end
function contract_center(a::Array{T, N}, b::Array{T, N}) where {T, N}
	(size(a) == size(b)) || throw(ArgumentError("a, b size mismatch"))
	a1 = reshape(a, length(a))
	b1 = reshape(b, length(b))
	@tensor r = a1[1] * b1[1]
	return r
end