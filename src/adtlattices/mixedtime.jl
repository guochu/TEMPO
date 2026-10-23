"""
    MixedADTLattice{O<:MixedFockOrdering} <: AbstractADTLattice{O}

ADT lattice on the mixed (Kadanoff-Baym, real-time + imaginary-time) contour: scalar type `ComplexF64`,
containing the `:+`, `:-` and `:τ` branches. The lattice discretizes the real-time interval `[0, t]` into `Nt` steps
(size `δt`) and the imaginary-time interval `[0, β]` into `Nτ` steps (size `δτ`).

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
(the extra one is the `τ=0` boundary point at lattice index 1, glued to the `:-` branch start by `boundarycondition!`;
the `τ=β` point at lattice index `kτ` is glued to the `:+` branch start); likewise each real-time branch stores
`kt = Nt + 1` points (the extra one is the `t=0` boundary point at lattice index 1, glued to the imaginary-time branch
ends). On the mixed contour, the correlation-function index `i` on **every** branch (`:τ`, `:+`, `:-`) maps to the
lattice index `i+1` — the same `lattice[BosonicContourIndex(i+1, branch)]` convention as the imaginary-time-only lattice: the lattice index 1
of each branch is the boundary/junction point and receives no influence-functional gates, while the lattice indices
`2..kτ` (`2..kt`) carry the correlation indices `1..Nτ` (`1..Nt`).

This convention is applied uniformly in `partialif`/`partialif_naive`/`hybriddynamics!` on the mixed contour
(see `src/influencefunctional/partialif/util.jl`, `partialif.jl` and `mixedtime.jl`).
"""
struct MixedADTLattice{O<:MixedFockOrdering} <: AbstractADTLattice{O}
	δt::Float64
	Nt::Int
	δτ::Float64
	Nτ::Int
	d::Int
	ordering::O

	MixedADTLattice(δt::Real, N::Int, δτ::Real, Ni::Int, d::Int, ordering::MixedFockOrdering) = new{typeof(ordering)}(
								convert(Float64, δt), N, convert(Float64, δτ), Ni, d, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
"""
    MixedADTLattice(; δt::Real, Nt::Int, δτ::Real, Nτ::Int, d::Int=2, ordering::MixedFockOrdering=M2M1_m1M1m2M2())

Convenience constructor for a mixed ADT lattice.
"""
MixedADTLattice(; δt::Real, Nt::Int, δτ::Real, Nτ::Int, d::Int=2, ordering::MixedFockOrdering=M2M1_m1M1m2M2()) = MixedADTLattice(
							δt, Nt, δτ, Nτ, d, ordering)
Base.similar(x::MixedADTLattice; δt::Real=x.δt, Nt::Int=x.Nt, δτ::Real=x.δτ, Nτ::Int=x.Nτ, d::Int=x.d, ordering::MixedFockOrdering=x.ordering) = MixedADTLattice(
			δt, Nt, δτ, Nτ, d, ordering)

TO.scalartype(::Type{<:MixedADTLattice}) = ComplexF64
"""
    branches(::Type{<:MixedADTLattice})

Tuple of branch symbols of the mixed ADT lattice, `(:+, :-, :τ)`.
"""
branches(::Type{<:MixedADTLattice}) = (:+, :-, :τ)

function Base.getproperty(x::MixedADTLattice, s::Symbol)
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

Base.length(x::MixedADTLattice) = 2 * x.kt + x.kτ

# acending order for real branch, descending order for imag time
"""
    getindex(x::MixedADTLattice{<:M2M1_m1M1m2M2}, a::BosonicContourIndex)

Map the contour index `a` (time step `a.j`, branch `:+`, `:-` or `:τ`) to the global lattice position.
"""
function Base.getindex(x::MixedADTLattice{<:M2M1_m1M1m2M2}, a::BosonicContourIndex)
	b = branch(a)
	@boundscheck begin
		(b in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
		if b == :τ
			(1 <= a.j <= x.kτ) || throw(BoundsError(1:x.kτ, a.j))
		else
			(1 <= a.j <= x.kt) || throw(BoundsError(1:x.kt, a.j))
		end
	end

	if b == :+
		2*(a.j-1)+2 + x.kτ
	elseif b == :-
		2*(a.j-1)+1 + x.kτ
	else
		x.kτ-a.j + 1
	end
end



# key is timestep, conj, branch, band
"""
    indexmappings(lattice::MixedADTLattice)

Return a `Dict{Tuple{Int, Symbol}, Int}` mapping `(time step, branch)` to the global index of the mixed ADT lattice.
"""
function indexmappings(lattice::MixedADTLattice)
	r = Dict{Tuple{Int, Symbol}, Int}()
	for i in 1:lattice.kτ
		f = :τ
		r[(i, f)] = lattice[BosonicContourIndex(i, f)]
	end
	for i in 1:lattice.kt
		for f in (:+, :-)
			r[(i, f)] = lattice[BosonicContourIndex(i, f)]
		end
	end
	return r
end

