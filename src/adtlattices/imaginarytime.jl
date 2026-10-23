"""
    ImagADTLattice{O<:ImagFockOrdering} <: AbstractADTLattice{O}

ADT lattice on the imaginary-time (Matsubara) contour: scalar type `Float64`, containing only the `:τ` branch.
The lattice discretizes the imaginary-time interval `[0, β]` into `N` steps of size `δτ = β/N`.

# Fields
- `δτ::Float64`: imaginary-time step size.
- `d::Int`: physical dimension.
- `N::Int`: number of imaginary-time discretization steps.
- `ordering::O`: Fock ordering (default [`M2M1`](@ref)).

# Derived attributes
- `β = N * δτ`: inverse temperature; `T = 1/β`.
- `τs = 0:δτ:β`: imaginary-time grid.
- `k = kτ = N + 1`: number of lattice indices (including boundary points).
"""
struct ImagADTLattice{O<:ImagFockOrdering} <: AbstractADTLattice{O}
	δτ::Float64
	d::Int
	N::Int
	ordering::O

	ImagADTLattice(δτ::Real, d::Int, N::Int, ordering::ImagFockOrdering) = new{typeof(ordering)}(float(δτ), d, N, ordering)
end

"""
    ImagADTLattice(; δτ::Real, N::Int, d::Int=2, ordering::ImagFockOrdering=M2M1())

Convenience constructor for an imaginary-time ADT lattice.
"""
ImagADTLattice(; δτ::Real, N::Int, d::Int=2, ordering::ImagFockOrdering=M2M1()) = ImagADTLattice(δτ, d, N, ordering)
Base.similar(x::ImagADTLattice; δτ::Real=x.δτ, d::Int=x.d, N::Int=x.N, ordering::ImagFockOrdering=x.ordering) = ImagADTLattice(δτ, d, N, ordering)

TO.scalartype(::Type{<:ImagADTLattice}) = Float64
"""
    branches(::Type{<:ImagADTLattice})

Tuple of branch symbols of the imaginary-time ADT lattice, `(:τ,)`.
"""
branches(::Type{<:ImagADTLattice}) = (:τ,)
TimeOrderingStyle(x::ImagADTLattice) = ImaginaryTimeOrderingStyle(x)

Base.length(x::ImagADTLattice) = x.k
function Base.getproperty(x::ImagADTLattice, s::Symbol)
	if s == :τs
		return 0:x.δτ:x.N*x.δτ
	elseif s == :β
		return x.N * x.δτ
	elseif s == :T
		return 1 / x.β
	elseif s == :Nτ
		return x.N
	elseif (s == :k) || ( s == :kτ)
		return x.N + 1
	else
		getfield(x, s)
	end
end


"""
    getindex(x::ImagADTLattice{<:M2M1}, a::BosonicContourIndex)

Map the contour index `a` (time step `a.j`, branch `:τ`) to the global lattice position.
"""
function Base.getindex(x::ImagADTLattice{<:M2M1}, a::BosonicContourIndex)
	@boundscheck begin
		(1 <= a.j <= x.k) || throw(BoundsError(1:x.k, a.j))
		(branch(a) == :τ) || throw(ArgumentError("branch must be :τ"))
	end
	return x.k-a.j+1
end



"""
    indexmappings(lattice::ImagADTLattice)

Return a `Dict{Tuple{Int, Symbol}, Int}` mapping `(time step, branch)` to the global index of the ADT lattice.
"""
function indexmappings(lattice::ImagADTLattice)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.k
		r[(i, :τ)] = lattice[BosonicContourIndex(i, :τ)]
	end
	return r
end
