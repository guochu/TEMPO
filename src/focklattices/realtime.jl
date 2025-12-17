abstract type RealFockLattice{O<:RealFockOrdering} <: AbstractFockLattice{O} end
TO.scalartype(::Type{<:RealFockLattice}) = ComplexF64
branches(::Type{<:RealFockLattice}) = (:+, :-)
TimeOrderingStyle(x::RealFockLattice) = RealTimeOrderingStyle(x)

# k is the number of discretization, nbands is the number of bands
# pos is the position within a band
# k+1 due to the grassmann number on the boundary for the final trace

"""
	struct RealFockLattice1Order <: RealFockLattice

First order splitting of the real-time contour
"""
struct RealFockLattice1Order{O<:RealFockOrdering} <: RealFockLattice{O}
	δt::Float64
	d::Int
	N::Int
	ordering::O

	RealFockLattice1Order(δt::Real, d::Int, N::Int, ordering::RealFockOrdering) = new{typeof(ordering)}(float(δt), d, N, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
RealFockLattice1Order(; δt::Real, N::Int, d::Int=2, ordering::RealFockOrdering=M1m1N1n1()) = RealFockLattice1Order(δt, d, N, ordering)
Base.similar(x::RealFockLattice1Order; δt::Real=x.δt, d::Int=x.d, N::Int=x.N, ordering::RealFockOrdering=x.ordering) = RealFockLattice1Order(δt, d, N, ordering)
# similargrassmannlattice(x::RealFockLattice1Order, δt::Real=x.δt, bands::Int=x.bands, N::Int=x.N, 
# 						ordering::RealGrassmannOrdering=similargrassmannordering(x.ordering)) = GrassmannLattice(contour=:real, δt=δt, N=N, bands=bands, ordering=ordering)


function Base.getproperty(x::RealFockLattice, s::Symbol)
	if s == :t
		return x.N * x.δt
	elseif s == :ts
		return 0:x.δt:x.N*x.δt
	elseif s == :Nt
		return x.N
	elseif (s == :k) || (s == :kt)
		return x.N+1
	else
		getfield(x, s)
	end
end
Base.length(x::RealFockLattice) = 2*x.k

function RealFockLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return RealFockLattice1Order(; kwargs...)
	else
		error("Second orderr RealFockLattice not implemented")
	end
end


function index(x::RealFockLattice{<:M2m2M1m1}, i::Int; branch::Symbol=:+)
	@boundscheck begin
		(1 <= i <= x.k) || throw(BoundsError(1:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	ifelse(branch == :+, TL-2i+1, TL-2i+2)
end


# key is timestep, conj, branch, band
function indexmappings(lattice::RealFockLattice)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.k
		for f in (:+, :-)
			r[(i, f)] = index(lattice, i, branch=f)
		end
	end
	return r
end

