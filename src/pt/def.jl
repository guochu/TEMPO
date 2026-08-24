
"""
	ProcessTensor{T<:Number, R<:Real}

One-dimensional dense tensor network (`Dense1DTN`) type representing a finite matrix product operator (MPO), storing a column of rank-4 site tensors.

# Fields
- `data::Vector{Array{T,4}}`: list of site tensors
- `s::Vector{Union{Missing, Vector{R}}}`: singular (Schmidt) vectors; `missing` when uninitialized
- `scaling::Ref{Float64}`: overall scaling factor

# Examples
```julia
julia> h = ProcessTensor(4, d=2)   # construct a ProcessTensor with 4 sites and physical dimension 2 (identity operator)
julia> h = randompt(4, D=8)        # construct a random ProcessTensor with 4 sites and bond dimension 8
```
"""
struct ProcessTensor{T<:Number,  R<:Real} <: Dense1DTN{T}
	data::Vector{Array{T, 4}}
	s::Vector{Union{Missing, Vector{R}}}
	scaling::Ref{Float64}
"""
	ProcessTensor{T, R}(data::AbstractVector, svectors::Vector, scaling::Ref{R}) where {T<:Number, R<:Number}

Inner constructor of `ProcessTensor`; only supports operators with strictly conserved quantum numbers.

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
function ProcessTensor{T, R}(data::AbstractVector, svectors::Vector, scaling::Ref{R}) where {T<:Number, R<:Number}
	(R == real(T)) || throw(ArgumentError("scalar type for singular vectors must be real"))
	(length(data)+1 == length(svectors)) || throw(DimensionMismatch("length of singular vectors must be length of site tensors+1"))
	_check_mpo_space(data)
	new{T, R}(convert(Vector{Array{T, 4}}, data), convert(Vector{Union{Missing, Vector{R}}}, svectors), scaling)
end
end

function ProcessTensor{T, R}(data::Vector, scaling::Ref{R}) where {T<:Number, R<:Number}
	(R == real(T)) || throw(ArgumentError("scalar type for singular vectors must be real"))
	_check_mpo_space(data)
	svectors = Vector{Union{Missing, Vector{R}}}(undef, length(data)+1)
	svectors[1] = ones(space_l(data[1]))
	svectors[end] = ones(space_r(data[end]))
	return ProcessTensor{T, R}(convert(Vector{Array{T, 4}}, data), svectors, scaling)
end

function ProcessTensor(data::AbstractVector{<:DenseMPOTensor{T}}, svectors::AbstractVector; scaling::Real=1) where {T <: Number}
	R = real(T)
	return ProcessTensor{T, R}(data, svectors, Ref(convert(R, scaling)))
end 
function ProcessTensor(data::AbstractVector{<:DenseMPOTensor{T}}; scaling::Real=1) where {T <: Number}
	R = real(T)
	return ProcessTensor{T, R}(data, Ref(convert(R, scaling)))
end


function ProcessTensor(::Type{T}, ds::AbstractVector{Int}) where {T <: Number}
	data = [reshape(_eye(T, d), 1, d, 1, d) for d in ds]
	return ProcessTensor(data, scaling=1)
end
ProcessTensor(ds::AbstractVector{Int}) = ProcessTensor(Float64, ds)
ProcessTensor(::Type{T}, L::Int; d::Int=2) where {T <: Number} = ProcessTensor(T, [d for i in 1:L])
ProcessTensor(L::Int; d::Int=2) = ProcessTensor(Float64, L, d=d)

Base.copy(psi::ProcessTensor) = ProcessTensor(copy(psi.data), copy(psi.s), scaling=scaling(psi))
function Base.copy!(a::ProcessTensor, b::ProcessTensor)
	a.data .= b.data
	a.s .= b.s
	setscaling!(a, scaling(b))
	return a
end
function Base.complex(psi::ProcessTensor)
	if scalartype(psi) <: Real
		data = [complex(item) for item in psi.data]
		return ProcessTensor(data, psi.s, scaling=scaling(psi))
	end
	return psi
end

svectors_uninitialized(psi::ProcessTensor) = any(ismissing, psi.s)
function unset_svectors!(psi::ProcessTensor)
	psi.s[2:end-1] .= missing
	return psi
end


function increase_bond!(psi::ProcessTensor, D::Int)
	if bond_dimension(psi) < D
		L = length(psi)
		for i in 1:L
			sl = (i == 1) ? 1 : max(D, size(psi[i], 1))
			sr = (i == L) ? 1 : max(D, size(psi[i], 3))
			m = zeros(scalartype(psi), sl, size(psi[i], 2), sr, size(psi[i], 4))
			m[1:size(psi[i], 1), :, 1:size(psi[i], 3), :] .= psi[i]
			psi[i] = m
		end
	end
	return psi
end

# attributes

function _check_mpo_space(mpotensors::Vector)
	L = length(mpotensors)
	for i in 1:L-1
		(space_r(mpotensors[i]) == space_l(mpotensors[i+1])) || throw(DimensionMismatch())
	end
	# for m in mpotensors
	# 	(size(m, 2) == size(m, 4)) || throw(ArgumentError("physical dimension mismatch"))
	# end
	# boundaries should be dimension 
	(space_l(mpotensors[1]) == 1) || throw(DimensionMismatch())
	(space_r(mpotensors[L]) == 1) || throw(DimensionMismatch())
	return true
end

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
randompt(::Type{T}, L::Int; d::Int=2, D::Int) where {T<:Number} = randompt(T, [d for i in 1:L], D=D)
randompt(L::Int; kwargs...) = randompt(Float64, L; kwargs...)



