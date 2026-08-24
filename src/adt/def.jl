"""
	ADT{T<:Number, R<:Real}

One-dimensional dense tensor network (`Dense1DTN`) type representing a finite matrix product state (MPS).

`ADT` stores a column of rank-3 site tensors and represents a one-dimensional quantum state with open boundary conditions (e.g., the discretized influence functional in the TEMPO algorithm). Site tensor conventions:

- Dimension 1: left auxiliary (bond) index, with dimension 1 at the leftmost site
- Dimension 2: physical index, with entry 1 corresponding to state |0⟩ and entry 2 to state |1⟩
- Dimension 3: right auxiliary (bond) index, with dimension 1 at the rightmost site

# Fields
- `data::Vector{Array{T,3}}`: list of site tensors
- `s::Vector{Union{Missing, Vector{R}}}`: singular (Schmidt) vectors; `missing` when uninitialized
- `scaling::Ref{Float64}`: overall scaling factor

# Examples
```julia
julia> psi = ADT(4, d=2)          # construct an ADT with 4 sites and physical dimension 2
julia> psi = randomadt(4, D=8)    # construct a random ADT with 4 sites and bond dimension 8
```
"""
struct ADT{T<:Number, R<:Real} <: Dense1DTN{T}
	data::Vector{Array{T, 3}}
	s::Vector{Union{Missing, Vector{R}}}
	scaling::Ref{Float64}

function ADT{T, R}(data::AbstractVector, svectors::Vector, scaling::Ref{R}) where {T<:Number, R<:Number}
	(R == real(T)) || throw(ArgumentError("scalar type for singular vectors must be real"))
	(length(data)+1 == length(svectors)) || throw(DimensionMismatch("length of singular vectors must be length of site tensors+1"))
	_check_mps_space(data)
	new{T, R}(convert(Vector{Array{T, 3}}, data), convert(Vector{Union{Missing, Vector{R}}}, svectors), scaling)
end
end

function ADT{T, R}(data::Vector, scaling::Ref{R}) where {T<:Number, R<:Number}
	(R == real(T)) || throw(ArgumentError("scalar type for singular vectors must be real"))
	_check_mps_space(data)
	svectors = Vector{Union{Missing, Vector{R}}}(undef, length(data)+1)
	svectors[1] = ones(space_l(data[1]))
	svectors[end] = ones(space_r(data[end]))
	return ADT{T, R}(convert(Vector{Array{T, 3}}, data), svectors, scaling)
end

function ADT(data::AbstractVector{<:DenseMPSTensor{T}}, svectors::AbstractVector; scaling::Real=1) where {T <: Number}
	R = real(T)
	return ADT{T, R}(data, svectors, Ref(convert(R, scaling)))
end 
function ADT(data::AbstractVector{<:DenseMPSTensor{T}}; scaling::Real=1) where {T <: Number}
	R = real(T)
	return ADT{T, R}(data, Ref(convert(R, scaling)))
end

# function ADT(::Type{T}, L::Int) where {T <: Number}
# 	v = zeros(T, 1, 2, 1)
# 	v[1,1,1] = 1
# 	data = [copy(v) for i in 1:L]
# 	return ADT(data, scaling=1)
# end
"""
	ADT(::Type{T}, ds::AbstractVector{Int}) where {T<:Number}

Construct an `ADT` with bond dimension 1 and all entries equal to one, whose physical dimensions are given by `ds`.

# Arguments
- `T`: element type (e.g. `Float64`, `ComplexF64`)
- `ds::AbstractVector{Int}`: physical dimension of each site

# Returns
An `ADT` whose entries are all `one(T)`.
"""
function ADT(::Type{T}, ds::AbstractVector{Int}) where {T <: Number}
	data = [ones(T, 1, d, 1) for d in ds]
	return ADT(data, scaling=1)
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
ADT(::Type{T}, L::Int; d::Int=2) where {T <: Number} = ADT(T, [d for i in 1:L])
"""
	ADT(L::Int; d::Int=2)

Construct an all-ones `ADT` with `L` sites of physical dimension `d` and element type `Float64`.
"""
ADT(L::Int; d::Int=2) = ADT(Float64, L, d=d)

Base.copy(psi::ADT) = ADT(copy(psi.data), copy(psi.s), scaling=scaling(psi))
function Base.copy!(a::ADT, b::ADT)
	a.data .= b.data
	a.s .= b.s
	setscaling!(a, scaling(b))
	return a
end
function Base.complex(psi::ADT)
	if scalartype(psi) <: Real
		data = [complex(item) for item in psi.data]
		return ADT(data, psi.s, scaling=scaling(psi))
	end
	return psi
end

svectors_uninitialized(psi::ADT) = any(ismissing, psi.s)
function unset_svectors!(psi::ADT)
	psi.s[2:end-1] .= missing
	return psi
end

function _check_mps_space(mpstensors::Vector)
	L = length(mpstensors)
	for i in 1:L-1
		(space_r(mpstensors[i]) == space_l(mpstensors[i+1])) || throw(DimensionMismatch())
	end
	(space_l(mpstensors[1]) == 1) || throw(DimensionMismatch("left boundary should be size 1"))
	(space_r(mpstensors[L]) == 1) || throw(DimensionMismatch("right boundary should be size 1"))
	return true
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
function randomadt(::Type{T}, ds::AbstractVector{Int}; D::Int) where {T <: Number}
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
randomadt(::Type{T}, L::Int; D::Int, d::Int=2) where {T<:Number} = randomadt(T, [d for i in 1:L], D=D)
"""
	randomadt(L::Int; kwargs...)

Generate a random `ADT` with `L` sites and element type `Float64`.
"""
randomadt(L::Int; kwargs...) = randomadt(Float64, L; kwargs...)

function increase_bond!(psi::ADT, D::Int)
	if bond_dimension(psi) < D
		L = length(psi)
		for i in 1:L
			sl = (i == 1) ? 1 : max(D, size(psi[i], 1))
			sr = (i == L) ? 1 : max(D, size(psi[i], 3))
			m = zeros(scalartype(psi), sl, size(psi[i], 2), sr)
			m[1:size(psi[i], 1), :, 1:size(psi[i], 3)] .= psi[i]
			psi[i] = m
		end
	end
	return psi
end


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


function easy_swap!(x::ADT, bond::Int; trunc::TruncationScheme=DefaultTruncation)
	x[bond], x.s[bond+1], x[bond+1] = _swap_gate(x.s[bond], x[bond], x.s[bond+1], x[bond+1], trunc=trunc)
	# x[bond], x[bond+1] = _swap_gate(x[bond], x[bond+1], trunc=trunc)
	return x
end

function naive_swap!(x::ADT, bond::Int; trunc::TruncationScheme=DefaultTruncation)
	x[bond], x[bond+1] = _swap_gate(x[bond], x[bond+1], trunc=trunc)
	# x[bond], x[bond+1] = _swap_gate(x[bond], x[bond+1], trunc=trunc)
	return x
end


function _swap_gate(m1, m2; trunc)
	@tensor twositemps[1,4;2,5] := m1[1,2,3] * m2[3,4,5]
	u, s, v, err = tsvd!(twositemps; trunc=trunc)
	# return u, permute(s * v, (1,2), (3,))
	return u * s, permute(v, (1,2), (3,))
end

function _swap_gate(svectorj1, m1, svectorj2, m2; trunc::TruncationScheme)
	@tensor twositemps[1,4;2,5] := m1[1,2,3] * m2[3,4,5]
	# println(space(svectorj1, 2), " ", space(m1, 1))
	@tensor twositemps1[-1 -2; -3 -4] := svectorj1[-1, 1] * twositemps[1, -2, -3, -4]
	u, s, v, err = tsvd!(twositemps1, trunc=trunc)
	@tensor u[-1 -2; -3] = twositemps[-1,-2,1,2] * conj(v[-3,1,2])
	return u, s, permute(v, (1,2), (3,))
end