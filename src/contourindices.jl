abstract type AbstractFockLatticeIndex end
# band(x::AbstractLatticeIndex) = x.band

"""
	ContourIndex <: AbstractFockLatticeIndex

Index on the Keldysh contour, consisting of a discrete time step `j` and a branch `branch`, used to index positions on a PTLattice.

Branch meanings:
- `:τ`: imaginary-time branch
- `:+`: forward real-time branch
- `:-`: backward real-time branch

Main constructor: `ContourIndex(j::Int; branch::Symbol=:τ)`.
"""
struct ContourIndex <: AbstractFockLatticeIndex
	j::Int
	branch::Symbol


function ContourIndex(j::Int, branch::Symbol)
	(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
	new(j, branch)
end

end
"""
	ContourIndex(j::Int; branch::Symbol=:τ)

Construct a contour index, with default branch `:τ` (imaginary time).

# Throws
- `ArgumentError`: `branch` is not one of `:τ`, `:+`, `:-`
"""
ContourIndex(j::Int; branch::Symbol=:τ) = ContourIndex(j, branch)
branch(x::ContourIndex) = x.branch

Base.:(==)(a::ContourIndex, b::ContourIndex) = (a.j == b.j) && (branch(a) == branch(b))

function Base.:<(a::ContourIndex, b::ContourIndex)
	(a == b) && return false
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

function Base.:(<=)(a::ContourIndex, b::ContourIndex)
	if a == b
		return true
	elseif a < b
		return true
	else
		return false
	end
end
Base.:(>)(a::ContourIndex, b::ContourIndex) = !(a <= b)
Base.:(>=)(a::ContourIndex, b::ContourIndex) = !(a < b)