"""
    RealADTLattice{O<:RealFockOrdering} <: AbstractADTLattice{O}

ADT lattice on the real-time (Keldysh) contour: scalar type `ComplexF64`, containing the `:+` and `:-` branches.
The lattice discretizes the real-time interval `[0, t]` into `N` steps of size `δt = t/N`.

# Fields
- `δt::Float64`: real-time step size.
- `d::Int`: physical dimension.
- `N::Int`: number of real-time discretization steps.
- `ordering::O`: Fock ordering (default [`M2m2M1m1`](@ref)).

# Derived attributes
- `t = N * δt`: total evolution time; `ts = 0:δt:t`: real-time grid.
- `Nt = N`; `k = kt = N + 1`: number of lattice indices (including boundary points).
"""
struct RealADTLattice{O<:RealFockOrdering} <: AbstractADTLattice{O}
	δt::Float64
	d::Int
	N::Int
	ordering::O

	RealADTLattice(δt::Real, d::Int, N::Int, ordering::RealFockOrdering) = new{typeof(ordering)}(float(δt), d, N, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
"""
    RealADTLattice(; δt::Real, N::Int, d::Int=2, ordering::RealFockOrdering=M2m2M1m1())

Convenience constructor for a real-time ADT lattice.
"""
RealADTLattice(; δt::Real, N::Int, d::Int=2, ordering::RealFockOrdering=M2m2M1m1()) = RealADTLattice(δt, d, N, ordering)
Base.similar(x::RealADTLattice; δt::Real=x.δt, d::Int=x.d, N::Int=x.N, ordering::RealFockOrdering=x.ordering) = RealADTLattice(δt, d, N, ordering)

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
    getindex(x::RealADTLattice{<:M2m2M1m1}, a::BosonicContourIndex)

Map the contour index `a` (time step `a.j`, branch `:+` or `:-`) to the global lattice position.
"""
function Base.getindex(x::RealADTLattice{<:M2m2M1m1}, a::BosonicContourIndex)
	@boundscheck begin
		(1 <= a.j <= x.k) || throw(BoundsError(1:x.k, a.j))
		(branch(a) in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	ifelse(branch(a) == :+, TL-2a.j+1, TL-2a.j+2)
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
			r[(i, f)] = lattice[BosonicContourIndex(i, f)]
		end
	end
	return r
end

