# store mpo tensor as a matrix of mpotensors instead of a single mpotensor
abstract type AbstractSparseMPOTensor end
const SparseMPOTensorElement{M<:AbstractMatrix, T<:Number} = Union{M, T}


storage(m::AbstractSparseMPOTensor) = m.Os

Base.size(m::AbstractSparseMPOTensor) = size(storage(m))
Base.size(m::AbstractSparseMPOTensor, i::Int) = size(storage(m), i)
phydim(m::AbstractSparseMPOTensor) = m.d

function Base.getindex(m::AbstractSparseMPOTensor, j::Int, k::Int) 
	T = scalartype(m)
	r = getindex(storage(m), j, k)
	d = phydim(m)
	if isa(r, T)
		if r == zero(T)
			return zeros(T, d, d)
		else
			return r * _eye(T, d)
		end
	else
		return r
	end
end 
Base.lastindex(m::AbstractSparseMPOTensor) = lastindex(storage(m))
Base.lastindex(m::AbstractSparseMPOTensor, i::Int) = lastindex(storage(m), i)

Base.keys(x::AbstractSparseMPOTensor) = Iterators.filter(a->contains(x, a[1],a[2]),Iterators.product(1:size(x, 1),1:size(x, 2)))
opkeys(x::AbstractSparseMPOTensor) = Iterators.filter(a-> !isscal(x,a[1],a[2]),keys(x))
scalkeys(x::AbstractSparseMPOTensor) = Iterators.filter(a-> isscal(x,a[1],a[2]),keys(x))


"""
	isid(x::MPOTensor; kwargs...)
	isid(x::MPSBondTensor; kwargs...)

Check if given MPOTensor or MPSBondTensor is identity 
"""
isid(x::Number; kwargs...) = true, x

function isid(x::AbstractMatrix; kwargs...)
	id = _eye(scalartype(x), size(x, 1))
	return _is_prop_util(x, id; kwargs...)
end

function _is_prop_util(x, a; atol::Real=1.0e-14) 
	scal = dot(a,x)/dot(a,a)
	diff = x-scal*a
	scal = (scal ≈ 0.0) ? 0.0 : scal #shouldn't be necessary (and I don't think it is)
	return norm(diff)<atol,scal
end
