"""
    RealADTLattice{O<:RealFockOrdering} <: AbstractADTLattice{O}

Abstract ADT lattice type on the real-time (Keldysh) contour: scalar type `ComplexF64`, containing the `:+` and `:-` branches.
"""
abstract type RealADTLattice{O<:RealFockOrdering} <: AbstractADTLattice{O} end
TO.scalartype(::Type{<:RealADTLattice}) = ComplexF64
"""
    branches(::Type{<:RealADTLattice})

Tuple of branch symbols of the real-time ADT lattice, `(:+, :-)`.
"""
branches(::Type{<:RealADTLattice}) = (:+, :-)
TimeOrderingStyle(x::RealADTLattice) = RealTimeOrderingStyle(x)

# k is the number of discretization, nbands is the number of bands
# pos is the position within a band
# k+1 due to the grassmann number on the boundary for the final trace

"""
	struct RealADTLattice1Order{O<:RealFockOrdering} <: RealADTLattice{O}

First-order ADT lattice on the real-time contour: discretizes the real-time interval `[0, t]` into `N` steps of size `δt = t/N`.

# Fields
- `δt::Float64`: real-time step size.
- `d::Int`: physical dimension.
- `N::Int`: number of real-time discretization steps.
- `ordering::O`: Fock ordering (default [`M2m2M1m1`](@ref)).

# Derived attributes
- `t = N * δt`: total evolution time; `ts = 0:δt:t`: real-time grid.
- `Nt = N`; `k = kt = N + 1`: number of lattice indices (including boundary points).
"""
struct RealADTLattice1Order{O<:RealFockOrdering} <: RealADTLattice{O}
	δt::Float64
	d::Int
	N::Int
	ordering::O

	RealADTLattice1Order(δt::Real, d::Int, N::Int, ordering::RealFockOrdering) = new{typeof(ordering)}(float(δt), d, N, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
"""
    RealADTLattice1Order(; δt::Real, N::Int, d::Int=2, ordering::RealFockOrdering=M2m2M1m1())

Convenience constructor for a first-order real-time ADT lattice.
"""
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

"""
    RealADTLattice(; order::Int=1, kwargs...)

Construct a real-time ADT lattice; `order` specifies the splitting order (currently only first order is supported, i.e. `order == 1`).
"""
function RealADTLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return RealADTLattice1Order(; kwargs...)
	else
		error("Second order RealADTLattice not implemented")
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
"""
    indexmappings(lattice::RealADTLattice)

Return a `Dict{Tuple{Int, Symbol}, Int}` mapping `(time step, branch)` to the global index of the ADT lattice.
"""
function indexmappings(lattice::RealADTLattice)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.k
		for f in (:+, :-)
			r[(i, f)] = index(lattice, i, branch=f)
		end
	end
	return r
end

