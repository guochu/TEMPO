"""
	ADT{T<:Number, R<:Real}

One-dimensional dense tensor network (`Dense1DTN`) type representing a finite matrix product state (MPS).

`ADT` stores a column of rank-3 site tensors and represents a one-dimensional quantum state with open boundary conditions (e.g., the discretized influence functional in the TEMPO algorithm). The storage payload is a FiniteMPSAlgorithms `CanonicalMPS` (field `.parent`): site tensors, Schmidt values and the per-site scaling are all carried by the payload, on which the FiniteMPSAlgorithms algorithms operate in place; `ADT` itself provides TEMPO's constructors and contour-specific operations (`swap!`, `permute!`, ...). `.data` delegates to the payload's site-tensor vector (same meaning as in old versions); `.s` / `.scaling` delegate to the payload's Schmidt values / scaling.

Site tensor conventions (payload `CanonicalMPS`):

- Dimension 1: left auxiliary (bond) index, with dimension 1 at the leftmost site
- Dimension 2: physical index, with entry 1 corresponding to state |0⟩ and entry 2 to state |1⟩
- Dimension 3: right auxiliary (bond) index, with dimension 1 at the rightmost site

# Examples
```julia
julia> psi = ADT(4, d=2)          # construct an ADT with 4 sites and physical dimension 2
julia> psi = randomadt(4, D=8)    # construct a random ADT with 4 sites and bond dimension 8
```
"""
struct ADT{T<:Number, R<:Real} <: Dense1DTN{T}
	parent::CanonicalMPS{T, R}
	ADT{T, R}(parent::CanonicalMPS{T, R}) where {T<:Number, R<:Real} = new{T, R}(parent)
end

ADT(parent::CanonicalMPS{T, R}) where {T<:Number, R<:Real} = ADT{T, R}(parent)

# `.parent` 即内层的 CanonicalMPS payload；`.data` 委托到 payload 的站点张量
# 向量（与旧版本语义一致）；`.s` / `.scaling` 委托到 payload
function Base.getproperty(psi::ADT, s::Symbol)
	s === :parent && return getfield(psi, :parent)
	s === :data && return getfield(psi, :parent).data
	s === :s && return getfield(psi, :parent).s
	s === :scaling && return getfield(psi, :parent).scaling
	throw(ArgumentError("ADT has no property $s"))
end
Base.propertynames(::ADT) = (:parent, :data, :s, :scaling)

function ADT(data::AbstractVector{<:DenseMPSTensor{T}}, svectors::AbstractVector; scaling::Real=1) where {T<:Number}
	R = real(T)
	mps = CanonicalMPS{T, R}(convert(Vector{Array{T, 3}}, data), svectors, convert(R, scaling))
	return ADT(mps)
end
function ADT(data::AbstractVector{<:DenseMPSTensor{T}}; scaling::Real=1) where {T<:Number}
	mps = CanonicalMPS{T, real(T)}(convert(Vector{Array{T, 3}}, data), convert(real(T), scaling))
	return ADT(mps)
end

"""
	ADT(::Type{T}, ds::AbstractVector{Int}) where {T<:Number}

Construct an `ADT` with bond dimension 1 and all entries equal to one, whose physical dimensions are given by `ds`.

# Arguments
- `T`: element type (e.g. `Float64`, `ComplexF64`)
- `ds::AbstractVector{Int}`: physical dimension of each site

# Returns
An `ADT` whose entries are all `one(T)`.
"""
function ADT(::Type{T}, ds::AbstractVector{Int}) where {T<:Number}
	return ADT(CanonicalMPS(T, ds))
end
"""
	ADT(ds::AbstractVector{Int})

Construct an all-ones `ADT` with element type `Float64` and physical dimensions given by `ds`.
"""
ADT(ds::AbstractVector{Int}) = ADT(Float64, ds)
"""
	ADT(::Type{T}, L::Int; d::Int=2) where {T<:Number}

Construct an all-ones `ADT` with `L` sites, each of physical dimension `d`.
"""
ADT(::Type{T}, L::Int; d::Int=2) where {T<:Number} = ADT(T, [d for _ in 1:L])
"""
	ADT(L::Int; d::Int=2)

Construct an all-ones `ADT` with `L` sites of physical dimension `d` and element type `Float64`.
"""
ADT(L::Int; d::Int=2) = ADT(Float64, L; d=d)

Base.copy(psi::ADT) = ADT(copy(psi.parent))
Base.copy!(a::ADT, b::ADT) = (copy!(a.parent, b.parent); a)
function Base.complex(psi::ADT)
	scalartype(psi) <: Real || return psi
	return ADT(complex(psi.parent))
end

svectors_uninitialized(psi::ADT) = svectors_uninitialized(psi.parent)
function unset_svectors!(psi::ADT)
	unset_svectors!(psi.parent)
	return psi
end

# initializers
"""
	randomadt(::Type{T}, ds::AbstractVector{Int}; D::Int) where {T<:Number}

Generate a randomly initialized `ADT` (MPS) with physical dimensions given by `ds` and bond dimension `D` at every bond.

# Arguments
- `T`: element type (e.g. `Float64`, `ComplexF64`)
- `ds::AbstractVector{Int}`: physical dimension of each site
- `D::Int`: bond dimension

# Returns
An `ADT` with random tensor entries.

# Examples
```julia
julia> psi = randomadt(ComplexF64, [2, 2, 2], D=16)
```
"""
function randomadt(::Type{T}, ds::AbstractVector{Int}; D::Int) where {T<:Number}
	L = length(ds)
	mpstensors = Vector{Array{T, 3}}(undef, L)
	mpstensors[1] = randn(T, 1,ds[1],D)
	mpstensors[end] = randn(T, D, ds[end], 1)
	for i in 2:L-1
		mpstensors[i] = randn(T, D, ds[i], D)
	end
	return ADT(mpstensors)
end
"""
	randomadt(ds::AbstractVector{Int}; kwargs...)

Generate a random `ADT` with element type `Float64`; remaining arguments as in `randomadt(::Type{T}, ds; D)`.
"""
randomadt(ds::AbstractVector{Int}; kwargs...) = randomadt(Float64, ds; kwargs...)
"""
	randomadt(::Type{T}, L::Int; D::Int, d::Int=2) where {T<:Number}

Generate a random `ADT` with `L` sites, physical dimension `d`, and bond dimension `D`.
"""
randomadt(::Type{T}, L::Int; D::Int, d::Int=2) where {T<:Number} = randomadt(T, [d for _ in 1:L], D=D)
"""
	randomadt(L::Int; kwargs...)

Generate a random `ADT` with `L` sites and element type `Float64`.
"""
randomadt(L::Int; kwargs...) = randomadt(Float64, L; kwargs...)

"""
	changebond!(psi::ADT, D::Int)

Bring the bond profile of the MPS to `min(D, feasible)`: bonds larger than the target
are shrunk by slicing the leading bond indices, smaller bonds are grown by zero padding
(the represented state is unchanged), followed by a no-truncation re-canonicalization.
Delegates to FiniteMPSAlgorithms' `changebond!(::CanonicalMPS; D)`; replaces the old
grow-only `increase_bond!`.
"""
changebond!(psi::ADT, D::Int) = (changebond!(psi.parent; D); psi)


# check is canonical
"""
	isleftcanonical(a::ADT; kwargs...)

Check whether all site tensors of the `ADT` are left-canonical.

# Returns
`true` if every site tensor satisfies the left-canonical condition (keyword arguments are passed to `isapprox` as tolerances).
"""
isleftcanonical(a::ADT; kwargs...) = all(x->isleftcanonical(x; kwargs...), a.data)
"""
	isrightcanonical(a::ADT; kwargs...)

Check whether all site tensors of the `ADT` are right-canonical.

# Returns
`true` if every site tensor satisfies the right-canonical condition (keyword arguments are passed to `isapprox` as tolerances).
"""
isrightcanonical(a::ADT; kwargs...) = all(x->isrightcanonical(x; kwargs...), a.data)

"""
	iscanonical(psi::ADT; kwargs...)

Check whether the MPS is in canonical form: all site tensors are right-canonical and the singular vectors are the correct Schmidt numbers.

# Returns
`true` if all of the following hold:
- all site tensors are right-canonical;
- the singular vectors are initialized (not cleared by `unset_svectors!`);
- for every bond, the squared singular vectors agree with the corresponding left contraction environment.

The canonical form improves the numerical stability of time evolution and enables efficient computation of observables for unitary systems.
"""
function iscanonical(psi::ADT; kwargs...)
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

# swap gate：委托给 FiniteMPSAlgorithms 的 CanonicalMPS `swap!`（Hastings
# 更新，需要时自动完成规范化；右规范形式与键谱在截断误差内保持）
function swap!(x::ADT, bond::Int; trunc::TruncationScheme=DefaultITruncation)
	swap!(x.parent, bond; trunc)
	return x
end
