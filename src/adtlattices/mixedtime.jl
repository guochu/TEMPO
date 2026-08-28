"""
    MixedADTLattice{O<:MixedFockOrdering} <: AbstractADTLattice{O}

Abstract ADT lattice type on the mixed (Kadanoff-Baym, real-time + imaginary-time) contour: scalar type `ComplexF64`,
containing the `:+`, `:-` and `:τ` branches.
"""
abstract type MixedADTLattice{O<:MixedFockOrdering} <: AbstractADTLattice{O} end
TO.scalartype(::Type{<:MixedADTLattice}) = ComplexF64
"""
    branches(::Type{<:MixedADTLattice})

Tuple of branch symbols of the mixed ADT lattice, `(:+, :-, :τ)`.
"""
branches(::Type{<:MixedADTLattice}) = (:+, :-, :τ)

"""
	struct MixedADTLattice1Order{O<:MixedFockOrdering} <: MixedADTLattice{O}

First-order ADT lattice on the mixed contour: discretizes the real-time interval `[0, t]` into `Nt` steps (size `δt`) and
the imaginary-time interval `[0, β]` into `Nτ` steps (size `δτ`).

# Fields
- `δt::Float64`: real-time step size.
- `Nt::Int`: number of real-time discretization steps.
- `δτ::Float64`: imaginary-time step size.
- `Nτ::Int`: number of imaginary-time discretization steps.
- `d::Int`: physical dimension.
- `ordering::O`: Fock ordering (default [`M2M1_m1M1m2M2`](@ref)).

# Derived attributes
- `t = Nt * δt`: total evolution time; `β = Nτ * δτ`: inverse temperature; `T = 1/β`.
- `ts = 0:δt:t`: real-time grid; `τs = 0:δτ:β`: imaginary-time grid.
- `kt = Nt + 1`, `kτ = Nτ + 1`: number of lattice indices per branch (including boundary points).

# Indexing convention

Although the imaginary-time evolution only has `Nτ` steps, the lattice stores `kτ = Nτ + 1` imaginary-time points
(the extra one is the boundary point, which the boundary condition connects to the real-time ends).
Consequently, every operation that maps a correlation-function index to a lattice position on the mixed contour
shifts the indices by **one time step** whenever the operation involves the imaginary-time branch:

- a correlation index `i` on the `:τ` branch maps to the lattice index `i+1`, i.e. to position `kτ - i` (the same
  `index(lattice, i+1)` convention as the imaginary-time-only lattice);
- when the partial-IF row lies on the `:τ` branch, **all** column indices `j` (including real-time-branch columns)
  are likewise taken at `j+1`;
- rows and columns on the real-time branches (`:+`, `:-`) use no offset (`index(lattice, i)`).

This convention is applied uniformly in `partialif`/`partialif_naive`/`hybriddynamics!` on the mixed contour
(see `src/influencefunctional/partialif/util.jl` and `src/influencefunctional/partialif/mixedtime.jl`).
"""
struct MixedADTLattice1Order{O<:MixedFockOrdering} <: MixedADTLattice{O}
	δt::Float64
	Nt::Int
	δτ::Float64
	Nτ::Int
	d::Int
	ordering::O

	MixedADTLattice1Order(δt::Real, N::Int, δτ::Real, Ni::Int, d::Int, ordering::MixedFockOrdering) = new{typeof(ordering)}(
								convert(Float64, δt), N, convert(Float64, δτ), Ni, d, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
"""
    MixedADTLattice1Order(; δt::Real, Nt::Int, δτ::Real, Nτ::Int, d::Int=2, ordering::MixedFockOrdering=M2M1_m1M1m2M2())

Convenience constructor for a first-order mixed ADT lattice.
"""
MixedADTLattice1Order(; δt::Real, Nt::Int, δτ::Real, Nτ::Int, d::Int=2, ordering::MixedFockOrdering=M2M1_m1M1m2M2()) = MixedADTLattice1Order(
							δt, Nt, δτ, Nτ, d, ordering)
Base.similar(x::MixedADTLattice1Order; δt::Real=x.δt, Nt::Int=x.Nt, δτ::Real=x.δτ, Nτ::Int=x.Nτ, d::Int=x.d, ordering::MixedFockOrdering=x.ordering) = MixedADTLattice1Order(
			δt, Nt, δτ, Nτ, d, ordering)


"""
    MixedADTLattice(; order::Int=1, kwargs...)

Construct a mixed ADT lattice; `order` specifies the splitting order (currently only first order is supported, i.e. `order == 1`).
"""
function MixedADTLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return MixedADTLattice1Order(; kwargs...)
	else
		error("Second orderr MixedGrassmannLattice not implemented")
	end
end

function Base.getproperty(x::MixedADTLattice1Order, s::Symbol)
	if s == :t
		return x.Nt * x.δt
	elseif s == :β
		return x.Nτ * x.δτ
	elseif s == :T
		return 1 / x.β
	elseif s == :ts
		return 0:x.δt:x.t
	elseif s == :τs
		return 0:x.δτ:x.β
	elseif s == :kt
		return x.Nt + 1
	elseif s == :kτ
		return x.Nτ + 1
	else
		getfield(x, s)
	end
end

Base.length(x::MixedADTLattice1Order) = 2 * x.kt + x.kτ

# acending order for real branch, descending order for imag time
function index(x::MixedADTLattice1Order{<:M2M1_m1M1m2M2}, i::Int; branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
		if branch == :τ
			(1 <= i <= x.kτ) || throw(BoundsError(1:x.kτ, i))
		else
			(1 <= i <= x.kt) || throw(BoundsError(1:x.kt, i))
		end
	end

	if branch == :+
		2*(i-1)+2 + x.kτ
	elseif branch == :-
		2*(i-1)+1 + x.kτ
	else
		x.kτ-i + 1
	end
end



# key is timestep, conj, branch, band
"""
    indexmappings(lattice::MixedADTLattice1Order)

Return a `Dict{Tuple{Int, Symbol}, Int}` mapping `(time step, branch)` to the global index of the mixed ADT lattice.
"""
function indexmappings(lattice::MixedADTLattice1Order)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.kτ
		f = :τ
		r[(i, f)] = index(lattice, i, branch=f)
	end
	for i in 1:lattice.kt
		for f in (:+, :-)
			r[(i, f)] = index(lattice, i, branch=f)
		end
	end
	return r
end

