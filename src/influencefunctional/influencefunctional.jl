abstract type InfluenceFunctionalAlgorithm end
"""
	PartialIF(; trunc::TruncateDimCutoff=Defaults.alg_trunc_if())

Construct the full influence functional as a product of partial influence functionals (partial IFs), each with bond dimension D=2.

See [SciPost Phys. Core 7, 063 (2024)].

# Arguments
- `trunc::TruncateDimCutoff=Defaults.alg_trunc_if()`: truncation scheme used when multiplying the partial IFs successively into the full IF.
"""
@kwdef struct PartialIF <: InfluenceFunctionalAlgorithm
	trunc::TruncateDimCutoff = Defaults.alg_trunc_if()
end

"""
	XTRGIF(; alg_expan, alg_evo, alg_mult, k, fast, verbosity)

Construct the influence functional as a translationally invariant MPO: first build the differential influence functional of width dt/2^k, then repeatedly square it k times (`fast=true`) via the tree bipartition scheme to obtain the full-length influence functional.

See [SciPost Phys. Core 7, 063 (2024)].

# Fields
- `alg_expan::ExponentialExpansionAlgorithm=Defaults.alg_expan()`: exponential (Prony) expansion algorithm for the bath correlation function.
- `alg_evo::TimeEvoMPOAlgorithm=Defaults.alg_evo()`: time-evolution algorithm for the differential influence functional (`FirstOrderStepper` or `ComplexStepper`).
- `alg_mult::MPSAlgorithm=Defaults.alg_mult()`: MPO multiplication (compression) algorithm (`DMRG1`, `DMRG2` or `SVDCompression`).
- `k::Int=5`: number of tree bipartition iterations, corresponding to 2^k time steps.
- `fast::Bool=true`: if `true`, use the fast tree bipartition scheme (k multiplications); otherwise use the sequential scheme (2^k-1 multiplications).
- `verbosity::Int=0`: verbosity level of the output.
"""
@kwdef struct XTRGIF{T<:ExponentialExpansionAlgorithm, E<:TimeEvoMPOAlgorithm, M<:MPSAlgorithm} <: InfluenceFunctionalAlgorithm
	alg_expan::T = Defaults.alg_expan()
	alg_evo::E = Defaults.alg_evo()
	alg_mult::M = Defaults.alg_mult()
	k::Int = 5
	fast::Bool = true
	verbosity::Int = 0
end

function Base.getproperty(x::XTRGIF, s::Symbol)
	if s == :trunc
		return x.alg_mult.trunc
	else
		getfield(x, s)
	end
end

"""
	AbstractTDVPIFAlgorithm <: InfluenceFunctionalAlgorithm

Abstract base of the TDVP-flow influence functional algorithms ([`TDVPIF`](@ref),
[`TDVP2IF`](@ref)): the influence functional is viewed as the "equilibrium state"
IF = exp(H) of the influence operator H (the MPO returned by `influenceoperators`),
and computed by evolving the identity influence functional (β = 0) along the
imaginary-time flow dz/dτ = H·z from τ = 0 to τ = 1.

On real-time lattices (`RealADTLattice` with `AdditiveHyb`, `RealPTLattice` with `GeneralHybStyle`) the influence operator driving the flow is the sum (direct sum) of the 4 branch MPOs returned by `influenceoperators`, since the site-wise product algebra satisfies e^a∘e^b = e^{a+b}.

Both algorithms are recommended for constructing the influence functional on pure imaginary-time or pure real-time contours. Mixed-contour lattices (`MixedADTLattice` / `MixedPTLattice`) are currently not supported; use `PartialIF` there.
"""
abstract type AbstractTDVPIFAlgorithm <: InfluenceFunctionalAlgorithm end

"""
	TDVPIF(; alg_expan, trunc, δ, verbosity)

Construct the influence functional with a second-order **single-site** TDVP imaginary-time flow: the influence functional is viewed as the "equilibrium state" IF = exp(H) of the influence operator H, and is computed by evolving the identity influence functional (β = 0) along the flow dz/dτ = H·z from τ = 0 to τ = 1.

Each flow step is one forward-backward TDVP sweep: the center tensor is evolved by +δτ/2 (Krylov exponentiation of the local effective map) and the bond matrix by -δτ/2, with a full +δτ step at the turning site. The state is zero-padded to bond dimension `trunc.D` before the flow, so the bond dimension grows with the correlations up to `trunc.D`.

# Fields
- `alg_expan::ExponentialExpansionAlgorithm=Defaults.alg_expan()`: exponential (Prony) expansion algorithm for the bath correlation function.
- `trunc::TruncateDimCutoff=Defaults.alg_trunc_if()`: bond dimension of the flow manifold / final influence functional.
- `δ::Float64=0.1`: imaginary-time step of the flow (0 < δ ≤ 1, adjusted so that 1/δ is an integer).
- `verbosity::Int=0`: verbosity level of the output.

See [`AbstractTDVPIFAlgorithm`](@ref) for the shared conventions (influence operator on real-time lattices, supported contours).
"""
@kwdef struct TDVPIF{T<:ExponentialExpansionAlgorithm} <: AbstractTDVPIFAlgorithm
	alg_expan::T = Defaults.alg_expan()
	trunc::TruncateDimCutoff = Defaults.alg_trunc_if()
	δ::Float64 = 0.1                      # imaginary-time step of the flow (0 < δ ≤ 1, adjusted so that 1/δ is an integer)
	verbosity::Int = 0
	function TDVPIF(alg_expan::T, trunc, δ, verbosity) where {T}
		(0 < δ <= 1) || throw(ArgumentError("δ must be a real number in (0, 1], got $δ"))
		# adjust δ so that 1/δ is an integer number of steps
		n = max(round(Int, 1 / δ), 1)
		new{T}(alg_expan, trunc, 1 / n, verbosity)
	end
end

"""
	TDVP2IF(; alg_expan, trunc, δ, verbosity)

Construct the influence functional with a second-order **two-site** TDVP imaginary-time flow — the same flow dz/dτ = H·z from τ = 0 to τ = 1 as [`TDVPIF`](@ref), but driven by FiniteMPSAlgorithms' `TDVP2` / `HadamardTDVP2` sweeps: the neighboring site pair is evolved jointly by +δτ/2 (Krylov exponentiation of the two-site effective generator), re-split by a *truncating* SVD under `trunc`, followed by the single-site backward half-step -δτ/2.

Unlike the single-site variant the bond dimension **adapts during the flow**: the bond-dimension-1 identity guess is grown by the pair updates themselves up to the cap `trunc.D` (no zero-padding preparation), and trimmed where the spectrum allows. `trunc` is therefore the truncation scheme of the pair re-splits, which both caps and shapes the final bond profile.

# Fields
- `alg_expan::ExponentialExpansionAlgorithm=Defaults.alg_expan()`: exponential (Prony) expansion algorithm for the bath correlation function.
- `trunc::TruncateDimCutoff=Defaults.alg_trunc_if()`: truncation scheme (carrying the bond cap `D`) of the two-site re-splits; the bond dimension grows during the flow up to `trunc.D`.
- `δ::Float64=0.1`: imaginary-time step of the flow (0 < δ ≤ 1, adjusted so that 1/δ is an integer).
- `verbosity::Int=0`: verbosity level of the output.

See [`AbstractTDVPIFAlgorithm`](@ref) for the shared conventions (influence operator on real-time lattices, supported contours).
"""
@kwdef struct TDVP2IF{T<:ExponentialExpansionAlgorithm} <: AbstractTDVPIFAlgorithm
	alg_expan::T = Defaults.alg_expan()
	trunc::TruncateDimCutoff = Defaults.alg_trunc_if()
	δ::Float64 = 0.1                      # imaginary-time step of the flow (0 < δ ≤ 1, adjusted so that 1/δ is an integer)
	verbosity::Int = 0
	function TDVP2IF(alg_expan::T, trunc, δ, verbosity) where {T}
		(0 < δ <= 1) || throw(ArgumentError("δ must be a real number in (0, 1], got $δ"))
		# adjust δ so that 1/δ is an integer number of steps
		n = max(round(Int, 1 / δ), 1)
		new{T}(alg_expan, trunc, 1 / n, verbosity)
	end
end




"""
	HybridizationStyle

Abstract type describing the form of the system-bath hybridization. The concrete styles are given by [`AdditiveHyb`](@ref) (additive/diagonal), [`NonAdditiveHyb`](@ref) (non-additive symmetric) and [`NonDiagonalHyb`](@ref) (off-diagonal).
"""
abstract type HybridizationStyle end


"""
	AdditiveHyb(op::AbstractVector{<:Real})
	AdditiveHyb(a::AbstractMatrix)

Additive (diagonal) system-bath hybridization style: the system-bath interaction takes the form V(a+a†), where the coupling operator V = Diagonal(op) is a diagonal matrix.

# Arguments
- `op::AbstractVector{<:Real}`: vector of diagonal elements, of length equal to the dimension of the system Hilbert space (i.e. `lattice.d`).
- `a::AbstractMatrix`: must be a diagonal matrix (otherwise an `ArgumentError` is thrown); its diagonal elements are extracted at construction.

# Notes
- `pairop(b::AdditiveHyb)` returns `(op, op)` with `op = Diagonal(b.op)`.
"""
struct AdditiveHyb <: HybridizationStyle
	op::Vector{Float64}
end

AdditiveHyb(x::AbstractVector{<:Real}) = AdditiveHyb(float(x))

phydim(b::AdditiveHyb) = length(b.op)

TO.scalartype(::Type{AdditiveHyb}) = Float64
"""
	pairop(hyb::HybridizationStyle)

Return the pair of system-bath coupling operators `(op1, op2)` such that the interaction term takes the form op1*a + op2*a′ (a is the fermionic annihilation operator).

- `AdditiveHyb`: returns `(Diagonal(op), Diagonal(op))`.
- `NonAdditiveHyb`: returns `(op, op)` (symmetric coupling).
- `NonDiagonalHyb`: returns `(op, op')` (op and its transpose couple to a and a′, respectively).
"""
function pairop(b::AdditiveHyb)
	op = Matrix(Diagonal(b.op))
	return op, op
end
function AdditiveHyb(a::AbstractMatrix)
	isdiag(a) || throw(ArgumentError("AdditiveHyb requires diagonal matrix"))
	adiag = [a[i, i] for i in 1:size(a, 1)]
	return AdditiveHyb(adiag)
end

abstract type GeneralHybStyle <: HybridizationStyle end

"""
	NonAdditiveHyb(op::AbstractMatrix{T}) where {T<:Number}

Non-additive (symmetric) system-bath hybridization style: the system-bath interaction takes the form op*(a+a'), i.e. the system operator op couples to both the fermionic annihilation operator a and the creation operator a'.

# Arguments
- `op::AbstractMatrix{T}`: square matrix (n x n) that must be Hermitian (otherwise an `ArgumentError` is thrown), where n equals the dimension of the system Hilbert space.

# Notes
- Can be converted from `AdditiveHyb` via `NonAdditiveHyb(hyb::AdditiveHyb)`.
- `pairop(b::NonAdditiveHyb)` returns `(op, op)`.
"""
struct NonAdditiveHyb{T<:Number} <: GeneralHybStyle
	op::Matrix{T}

function NonAdditiveHyb{T}(op::AbstractMatrix) where {T<:Number}
	(size(op, 1) == size(op, 2)) || throw(ArgumentError("square matrix expected"))
	ishermitian(op) || throw(ArgumentError("Hermitian matrix required"))
	new{T}(convert(Matrix{T}, op))
end
end
NonAdditiveHyb(a::AbstractMatrix{T}) where {T<:Number} = NonAdditiveHyb{T}(a)
NonAdditiveHyb(hyb::AdditiveHyb) = NonAdditiveHyb(Matrix(Diagonal(hyb.op)))

phydim(b::NonAdditiveHyb) = size(b.op, 1)

TO.scalartype(::Type{NonAdditiveHyb{T}}) where T = T
pairop(b::NonAdditiveHyb) = b.op, b.op

"""
	NonDiagonalHyb(op::AbstractMatrix{T}) where {T<:Number}

Off-diagonal system-bath hybridization style: the system-bath interaction takes the form op*a + op'*a', i.e. the system operator op couples to a and its transpose couples to a'.

# Arguments
- `op::AbstractMatrix{T}`: square matrix (n x n), where n equals the dimension of the system Hilbert space. Unlike `NonAdditiveHyb`, no Hermiticity is required here, so asymmetric (e.g. chiral) couplings can be described.

# Notes
- Can be converted from other styles via `NonDiagonalHyb(hyb::AdditiveHyb)` or `NonDiagonalHyb(hyb::NonAdditiveHyb)`.
- `pairop(b::NonDiagonalHyb)` returns `(op, op')`.
"""
struct NonDiagonalHyb{T<:Number} <: GeneralHybStyle
	op::Matrix{T}

function NonDiagonalHyb{T}(op::AbstractMatrix) where {T<:Number}
	(size(op, 1) == size(op, 2)) || throw(ArgumentError("square matrix expected"))
	new{T}(convert(Matrix{T}, op))
end
end
NonDiagonalHyb(a::AbstractMatrix{T}) where {T<:Number} = NonDiagonalHyb{T}(a)
NonDiagonalHyb(hyb::AdditiveHyb) = NonDiagonalHyb(Matrix(Diagonal(hyb.op)))
NonDiagonalHyb(hyb::NonAdditiveHyb) = NonDiagonalHyb(hyb.op)

phydim(b::NonDiagonalHyb) = size(b.op, 1)

TO.scalartype(::Type{NonDiagonalHyb{T}}) where T = T
pairop(b::NonDiagonalHyb) = b.op, b.op'

include("partialif/partialif.jl")
include("ptpartialif/ptpartialif.jl")
include("ttiif/ttiif.jl")
include("tdvpif/tdvpif.jl")
