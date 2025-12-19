abstract type AbstractFockLatticeIndex end
# band(x::AbstractLatticeIndex) = x.band

struct ContourIndex <: AbstractFockLatticeIndex
	j::Int
	branch::Symbol


function ContourIndex(j::Int, branch::Symbol)
	(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
	new(j, branch)
end

end
ContourIndex(j::Int; branch::Symbol=:τ) = ContourIndex(j, branch)
branch(x::ContourIndex) = x.branch

Base.:(==)(a::ContourIndex, b::ContourIndex) = (a.j == b.j) && (branch(a) == branch(b))

function Base.:<(a::ContourIndex, b::ContourIndex)
	(a == b) && throw(ArgumentError("the two indices are the same"))
	if branch(a) == :+
		if branch(b) == :+
			return a.j < b.j
		else
			return true
		end
	elseif branch(a) == :-
		if branch(b) == :+
			return false
		elseif branch(b) == :-
			return a.j > b.j
		else
			return true
		end
	else
		if branch(b) == :τ
			return a.j < b.j
		else
			return false
		end
	end
end


struct ContourOperator{T<:Number}
	indices::Vector{ContourIndex}
	ops::Vector{Matrix{T}}
end

function ContourOperator(idx::AbstractVector{ContourIndex}, data::AbstractVector{<:AbstractMatrix}) 
	(length(Set(idx)) == length(idx)) || throw(ArgumentError("multiple op with the same ContourIndex not allowed"))
	(length(idx) == length(data)) || throw(DimensionMismatch("op size mismatch with number of ContourIndexes"))
	return ContourOperator(idx, data)
end
ContourOperator(p::ContourIndex, data::AbstractMatrix) = ContourOperator([p], [data])

TO.scalartype(::Type{ContourOperator{T}}) where {T} = T

