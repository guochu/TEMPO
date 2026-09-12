abstract type AbstractImpurityOperator end

"""
    AbstractImpurityHamiltonian <: AbstractImpurityOperator

Supertype of the (unitary) impurity Hamiltonian models: constant
([`ImpurityHamiltonian`](@ref)), quenched ([`QuenchedImpurityHamiltonian`](@ref))
and time-dependent ([`TdImpurityHamiltonian`](@ref)). Interface: `propagator`
(constant per branch for the first two; per-step for the latter), `phydim`.
"""
abstract type AbstractImpurityHamiltonian <: AbstractImpurityOperator end


"""
    ImpurityHamiltonian{M<:AbstractMatrix} <: AbstractImpurityHamiltonian

Hamiltonian model of the impurity system, wrapping a `d × d` matrix `m`, where `d` is the physical dimension.
"""
struct ImpurityHamiltonian{M<:AbstractMatrix} <: AbstractImpurityHamiltonian
	m::M
end
propagator(h::ImpurityHamiltonian, lat, b::Symbol) = _get_propagator(h.m, lat, b)
propagator(h::ImpurityHamiltonian, lat; branch::Symbol=:τ) = propagator(h, lat, branch)
phydim(h::ImpurityHamiltonian) = size(h.m, 1)
"""
    ImpurityHamiltonian(d::Int)

Construct a zero-matrix Hamiltonian of physical dimension `d`.
"""
ImpurityHamiltonian(d::Int) = ImpurityHamiltonian(zeros(d, d))

TO.scalartype(::Type{ImpurityHamiltonian{M}}) where {M} = scalartype(M)


"""
    QuenchedImpurityHamiltonian{M1<:AbstractMatrix, M2<:AbstractMatrix} <: AbstractImpurityHamiltonian

Impurity Hamiltonian with a quench protocol: the impurity evolves with the
pre-quench Hamiltonian `hτ` on the imaginary-time branch (`:τ`) and with the
post-quench Hamiltonian `ht` on the real-time branches (`:+`/`:-`).
"""
struct QuenchedImpurityHamiltonian{M1<:AbstractMatrix, M2<:AbstractMatrix} <: AbstractImpurityHamiltonian
	hτ::M1
	ht::M2

	function QuenchedImpurityHamiltonian(hτ::AbstractMatrix, ht::AbstractMatrix)
		(size(hτ) == size(ht)) || throw(DimensionMismatch("hτ and ht must have the same size, got $(size(hτ)) and $(size(ht))"))
		return new{typeof(hτ), typeof(ht)}(hτ, ht)
	end
end
QuenchedImpurityHamiltonian(hτ::AbstractMatrix, ht::AbstractMatrix, d::Int) = begin
	(size(hτ) == (d, d) && size(ht) == (d, d)) || throw(DimensionMismatch("expected $(d)×$(d) matrices"))
	return QuenchedImpurityHamiltonian(hτ, ht)
end
propagator(h::QuenchedImpurityHamiltonian, lat, b::Symbol) = _get_propagator(b == :τ ? h.hτ : h.ht, lat, b)
propagator(h::QuenchedImpurityHamiltonian, lat; branch::Symbol=:τ) = propagator(h, lat, branch)
phydim(h::QuenchedImpurityHamiltonian) = size(h.hτ, 1)
TO.scalartype(::Type{QuenchedImpurityHamiltonian{M1, M2}}) where {M1, M2} = promote_type(scalartype(M1), scalartype(M2))


"""
    TdImpurityOp{M<:AbstractMatrix}

Time-dependent impurity Hamiltonian contribution: the matrix `m` multiplied
by the scalar `f(t)` gives the Hamiltonian contribution of this term at time `t`.
"""
struct TdImpurityOp{M<:AbstractMatrix}
	m::M
	f::Function
end
(x::TdImpurityOp)(t::Real) = x.f(t) .* x.m
TO.scalartype(::Type{TdImpurityOp{M}}) where {M} = scalartype(M)


"""
    TdImpurityHamiltonian{M1<:AbstractMatrix, M2<:AbstractMatrix} <: AbstractImpurityHamiltonian

Time-dependent impurity Hamiltonian. `hτ` drives the imaginary-time branch
(`:τ`), while the real-time branches (`:+`/`:-`) evolve with the
time-dependent Hamiltonian `ht + Σₖ (htt[k])(t)`: calling the model with a
time `t` returns the impurity Hamiltonian matrix of the real-time branches
at `t`. The propagator on the real-time branches is evaluated per step, with
`t` the left endpoint of the physical time interval of the step (forward
branch: `(j-1)δt`; backward branch: `(Nt-j)δt` for step `j` of `Nt`).
"""
struct TdImpurityHamiltonian{M1<:AbstractMatrix, M2<:AbstractMatrix, V<:AbstractVector{<:TdImpurityOp}} <: AbstractImpurityHamiltonian
	hτ::M1
	ht::M2
	htt::V

	function TdImpurityHamiltonian(hτ::AbstractMatrix, ht::AbstractMatrix, htt::AbstractVector{<:TdImpurityOp})
		(size(hτ) == size(ht)) || throw(DimensionMismatch("hτ and ht must have the same size, got $(size(hτ)) and $(size(ht))"))
		for op in htt
			(size(op.m) == size(hτ)) || throw(DimensionMismatch("TdImpurityOp matrix size $(size(op.m)) does not match the model size $(size(hτ))"))
		end
		return new{typeof(hτ), typeof(ht), typeof(htt)}(hτ, ht, htt)
	end
end
TdImpurityHamiltonian(hτ::AbstractMatrix, ht::AbstractMatrix, htt::AbstractVector{<:TdImpurityOp}, d::Int) = begin
	(size(hτ) == (d, d) && size(ht) == (d, d)) || throw(DimensionMismatch("expected $(d)×$(d) matrices"))
	return TdImpurityHamiltonian(hτ, ht, htt)
end
# the impurity Hamiltonian matrix of the real-time branches at time t
function (h::TdImpurityHamiltonian)(t::Real)
	H = copy(h.ht)
	for op in h.htt
		H .+= op(t)
	end
	return H
end
phydim(h::TdImpurityHamiltonian) = size(h.hτ, 1)
TO.scalartype(::Type{TdImpurityHamiltonian{M1, M2, V}}) where {M1, M2, V} = promote_type(scalartype(M1), scalartype(M2), scalartype(V))
TO.scalartype(::Type{V}) where {V<:AbstractVector{<:TdImpurityOp}} = scalartype(eltype(V))

# the imaginary-time branch evolves with the constant hτ; the real-time
# branches are time-dependent and require the step index
propagator(h::TdImpurityHamiltonian, lat, b::Symbol) = begin
	(b == :τ) || throw(ArgumentError("the propagator of a time-dependent model on the real-time branches is step dependent: use sysdynamics (which evaluates the per-step propagator)"))
	_get_propagator(h.hτ, lat, b)
end
propagator(h::TdImpurityHamiltonian, lat; branch::Symbol=:τ) = propagator(h, lat, branch)

"""
    propagator(model, lattice, branch, j, N)

Propagator of the model acting at step `j` of `N` steps on the given contour
branch. Defaults to the (time-independent) branch propagator; the
time-dependent model evaluates its Hamiltonian at the step time `t`
(forward branch: `t = (j-1)δt`; backward branch: `t = (N-j)δt`).
"""
propagator(h::AbstractImpurityHamiltonian, lat, b::Symbol, j::Int, N::Int) = propagator(h, lat, b)
function propagator(h::TdImpurityHamiltonian, lat, b::Symbol, j::Int, N::Int)
	(b in (:+, :-)) || return propagator(h, lat, b) # the :τ branch is time-independent
	t = (b == :+) ? (j - 1) * lat.δt : (N - j) * lat.δt
	return _get_propagator(h(t), lat, b)
end



# dissipative impurity
"""
    ImpurityLindbladian <: AbstractImpurityOperator

Lindblad-type dissipative model of the impurity system, wrapping a `d × d × d × d` tensor `m`
(in superoperator form, for dissipative real-time evolution).
"""
struct ImpurityLindbladian <: AbstractImpurityOperator
	m::Array{ComplexF64, 4}
end

phydim(h::ImpurityLindbladian) = size(h.m, 1)

"""
    ImpurityLindbladian(d::Int)

Construct a zero Lindblad superoperator of physical dimension `d`.
"""
ImpurityLindbladian(d::Int) = ImpurityLindbladian(zeros(ComplexF64, d, d, d, d))
TO.scalartype(::Type{ImpurityLindbladian}) = ComplexF64
"""
    ImpurityLindbladian(L::LindbladOperator)

Construct an impurity Lindblad model from a `LindbladOperator`.
"""
ImpurityLindbladian(L::LindbladOperator) = ImpurityLindbladian(L.m)
"""
    ImpurityLindbladian(H::AbstractMatrix, jumpops::Vector{<:AbstractMatrix})

Construct a Lindblad superoperator from the Hamiltonian `H` and the jump operators `jumpops`.
"""
ImpurityLindbladian(H::AbstractMatrix, jumpops::Vector{<:AbstractMatrix}) = ImpurityLindbladian(lindbladoperator(H, jumpops))
"""
    ImpurityLindbladian(h::ImpurityHamiltonian)

Construct the corresponding Lindblad model from an [`ImpurityHamiltonian`](@ref) without jump operators (purely unitary evolution).
"""
ImpurityLindbladian(h::ImpurityHamiltonian) = ImpurityLindbladian(lindbladoperator(h.m, []))