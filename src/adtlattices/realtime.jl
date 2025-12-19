abstract type RealADTLattice{O<:RealFockOrdering} <: AbstractADTLattice{O} end
TO.scalartype(::Type{<:RealADTLattice}) = ComplexF64
branches(::Type{<:RealADTLattice}) = (:+, :-)
TimeOrderingStyle(x::RealADTLattice) = RealTimeOrderingStyle(x)

# k is the number of discretization, nbands is the number of bands
# pos is the position within a band
# k+1 due to the grassmann number on the boundary for the final trace

"""
	struct RealADTLattice1Order <: RealADTLattice

First order splitting of the real-time contour
"""
struct RealADTLattice1Order{O<:RealFockOrdering} <: RealADTLattice{O}
	δt::Float64
	d::Int
	N::Int
	ordering::O

	RealADTLattice1Order(δt::Real, d::Int, N::Int, ordering::RealFockOrdering) = new{typeof(ordering)}(float(δt), d, N, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
RealADTLattice1Order(; δt::Real, N::Int, d::Int=2, ordering::RealFockOrdering=M2m2M1m1()) = RealADTLattice1Order(δt, d, N, ordering)
Base.similar(x::RealADTLattice1Order; δt::Real=x.δt, d::Int=x.d, N::Int=x.N, ordering::RealFockOrdering=x.ordering) = RealADTLattice1Order(δt, d, N, ordering)
# similargrassmannlattice(x::RealADTLattice1Order, δt::Real=x.δt, bands::Int=x.bands, N::Int=x.N, 
# 						ordering::RealGrassmannOrdering=similargrassmannordering(x.ordering)) = GrassmannLattice(contour=:real, δt=δt, N=N, bands=bands, ordering=ordering)


function Base.getproperty(x::RealADTLattice, s::Symbol)
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
Base.length(x::RealADTLattice) = 2*x.k

function RealADTLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return RealADTLattice1Order(; kwargs...)
	else
		error("Second orderr RealADTLattice not implemented")
	end
end


function index(x::RealADTLattice{<:M2m2M1m1}, i::Int; branch::Symbol=:+)
	@boundscheck begin
		(1 <= i <= x.k) || throw(BoundsError(1:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	ifelse(branch == :+, TL-2i+1, TL-2i+2)
end


# key is timestep, conj, branch, band
function indexmappings(lattice::RealADTLattice)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.k
		for f in (:+, :-)
			r[(i, f)] = index(lattice, i, branch=f)
		end
	end
	return r
end

