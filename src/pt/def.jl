
"""
	ProcessTensor{T<:Number, R<:Real}

One-dimensional dense tensor network (`Dense1DTN`) type representing a finite matrix product operator (MPO), storing a column of rank-4 site tensors.

The storage payload is a FiniteMPSAlgorithms `CanonicalMPO` (field `.parent`): site tensors, Schmidt values and the per-site scaling are all carried by the payload, on which the FiniteMPSAlgorithms algorithms operate in place; `ProcessTensor` itself provides TEMPO's constructors and contour-specific operations (`swap!`, `permute!`, ...). `.data` delegates to the payload's site-tensor vector (same meaning as in old versions); `.s` / `.scaling` delegate to the payload's Schmidt values / scaling.

# Fields
- `parent::CanonicalMPO`: the payload chain (site tensors `Vector{Array{T,4}}`)

# Examples
```julia
julia> h = ProcessTensor(4, d=2)   # construct a ProcessTensor with 4 sites and physical dimension 2 (identity operator)
julia> h = randompt(4, D=8)        # construct a random ProcessTensor with 4 sites and bond dimension 8
```
"""
struct ProcessTensor{T<:Number, R<:Real} <: Dense1DTN{T}
	parent::CanonicalMPO{T, R}
"""
	ProcessTensor(parent::CanonicalMPO)

Inner constructor wrapping a FiniteMPSAlgorithms `CanonicalMPO` payload.

Site tensor index convention (i denotes the input arrow, o the output arrow):

	    o
	    |
	    2
	o-1   3-i
	    4
	    |
	    i

Both left and right boundaries are vacuum (dimension 1). A non-vacuum right boundary corresponds to operators that do not conserve quantum numbers (e.g., a†); such operators must be represented in other MPO forms.
"""
	function ProcessTensor{T, R}(parent::CanonicalMPO{T, R}) where {T<:Number, R<:Real}
		new{T, R}(parent)
	end
end

ProcessTensor(parent::CanonicalMPO{T, R}) where {T<:Number, R<:Real} = ProcessTensor{T, R}(parent)

# `.parent` 即内层的 CanonicalMPO payload；`.data` 委托到 payload 的站点张量
# 向量（与旧版本语义一致）；`.s` / `.scaling` 委托到 payload
function Base.getproperty(psi::ProcessTensor, s::Symbol)
	s === :parent && return getfield(psi, :parent)
	s === :data && return getfield(psi, :parent).data
	s === :s && return getfield(psi, :parent).s
	s === :scaling && return getfield(psi, :parent).scaling
	throw(ArgumentError("ProcessTensor has no property $s"))
end
Base.propertynames(::ProcessTensor) = (:parent, :data, :s, :scaling)

function ProcessTensor(data::AbstractVector{<:DenseMPOTensor{T}}, svectors::AbstractVector; scaling::Real=1) where {T<:Number}
	R = real(T)
	mps = CanonicalMPO{T, R}(convert(Vector{Array{T, 4}}, data), svectors, convert(R, scaling))
	return ProcessTensor(mps)
end
function ProcessTensor(data::AbstractVector{<:DenseMPOTensor{T}}; scaling::Real=1) where {T<:Number}
	mps = CanonicalMPO{T, real(T)}(convert(Vector{Array{T, 4}}, data), convert(real(T), scaling))
	return ProcessTensor(mps)
end


function ProcessTensor(::Type{T}, ds::AbstractVector{Int}) where {T<:Number}
	return ProcessTensor(CanonicalMPO(T, ds))
end
ProcessTensor(ds::AbstractVector{Int}) = ProcessTensor(Float64, ds)
ProcessTensor(::Type{T}, L::Int; d::Int=2) where {T<:Number} = ProcessTensor(T, [d for _ in 1:L])
ProcessTensor(L::Int; d::Int=2) = ProcessTensor(Float64, L; d=d)

Base.copy(psi::ProcessTensor) = ProcessTensor(copy(psi.parent))
Base.copy!(a::ProcessTensor, b::ProcessTensor) = (copy!(a.parent, b.parent); a)
function Base.complex(psi::ProcessTensor)
	scalartype(psi) <: Real || return psi
	return ProcessTensor(complex(psi.parent))
end

svectors_uninitialized(psi::ProcessTensor) = svectors_uninitialized(psi.parent)
function unset_svectors!(psi::ProcessTensor)
	unset_svectors!(psi.parent)
	return psi
end


"""
	changebond!(g::ProcessTensor, D::Int)

Bring the bond profile of the MPO to `min(D, feasible)`: bonds larger than the target
are shrunk by slicing the leading bond indices, smaller bonds are grown by zero padding
(the represented operator is unchanged; the chain is never re-gauged in place).
Delegates to FiniteMPSAlgorithms' `changebond!(::AbstractMPO; D)`; replaces the old
grow-only `increase_bond!`.
"""
changebond!(g::ProcessTensor, D::Int) = (changebond!(g.parent; D); g)

"""
	isleftcanonical(a::ProcessTensor; kwargs...)

Check whether all site tensors of the `ProcessTensor` are left-canonical.

# Returns
`true` if every site tensor satisfies the left-canonical condition (keyword arguments are passed to `isapprox` as tolerances).
"""
isleftcanonical(a::ProcessTensor; kwargs...) = all(x->isleftcanonical(x; kwargs...), a.data)
"""
	isrightcanonical(a::ProcessTensor; kwargs...)

Check whether all site tensors of the `ProcessTensor` are right-canonical.

# Returns
`true` if every site tensor satisfies the right-canonical condition (keyword arguments are passed to `isapprox` as tolerances).
"""
isrightcanonical(a::ProcessTensor; kwargs...) = all(x->isrightcanonical(x; kwargs...), a.data)

"""
	iscanonical(psi::ProcessTensor; kwargs...)

Check whether the MPO is in canonical form: all site tensors are right-canonical and the singular vectors are the correct Schmidt numbers.

# Returns
`true` if all of the following hold:
- all site tensors are right-canonical;
- the singular vectors are initialized (not cleared by `unset_svectors!`);
- for every bond, the squared singular vectors agree with the corresponding left contraction environment.

The canonical form improves the numerical stability of time evolution and enables efficient computation of observables for unitary systems.
"""
function iscanonical(psi::ProcessTensor; kwargs...)
	isrightcanonical(psi) || return false
	# we also check whether the singular vectors are the correct Schmidt numbers
	svectors_uninitialized(psi) && return false
	hold = l_LL(psi, psi)
	for i in 1:length(psi)-1
		hold = updateleft(hold, psi[i], psi[i])
		tmp = psi.s[i+1]
		isapprox(hold, Diagonal(tmp.^2); kwargs...) || return false
	end
	return true
end

# initializers

"""
	randompt(::Type{T}, ds::Vector{Int}; D::Int) where {T<:Number}

Generate a randomly initialized `ProcessTensor` (MPO) with physical dimensions given by `ds` and bond dimension `D` at every bond.

# Arguments
- `T`: element type (e.g. `Float64`, `ComplexF64`)
- `ds::Vector{Int}`: physical dimension of each site
- `D::Int`: bond dimension

# Returns
A `ProcessTensor` with random tensor entries.

# Examples
```julia
julia> h = randompt(ComplexF64, [2, 2, 2], D=16)
```
"""
function randompt(::Type{T}, ds::Vector{Int}; D::Int) where {T<:Number}
	L = length(ds)
	r = Vector{Array{T, 4}}(undef, L)
	r[1] = randn(T, 1, ds[1], D, ds[1])
	r[L] = randn(T, D, ds[L], 1, ds[L])
	for i in 2:L-1
		r[i] = randn(T, D, ds[i], D, ds[i])
	end
	return ProcessTensor(r)
end
randompt(::Type{T}, L::Int; d::Int=2, D::Int) where {T<:Number} = randompt(T, [d for _ in 1:L], D=D)
randompt(L::Int; kwargs...) = randompt(Float64, L; kwargs...)
