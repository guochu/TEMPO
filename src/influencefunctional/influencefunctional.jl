abstract type InfluenceFunctionalAlgorithm end
"""
	PartialIF(; trunc::TruncationDimCutoff=DefaultITruncation)

Construct the full influence functional as a product of partial influence functionals (partial IFs), each with bond dimension D=2.

See [SciPost Phys. Core 7, 063 (2024)].

# Arguments
- `trunc::TruncationDimCutoff=DefaultITruncation`: truncation scheme used when multiplying the partial IFs successively into the full IF.
"""
struct PartialIF <: InfluenceFunctionalAlgorithm 
	trunc::TruncationDimCutoff
end
"""
	PartialIF(; trunc::TruncationDimCutoff=DefaultITruncation)

Keyword constructor for `PartialIF`; the truncation scheme is specified by `trunc`.
"""
PartialIF(; trunc::TruncationDimCutoff=DefaultITruncation) = PartialIF(trunc)

"""
	TranslationInvariantIF(; algexpan, algevo, algmult, k, fast, verbosity)

Construct the influence functional as a translationally invariant MPO: first build the differential influence functional of width dt/2^k, then repeatedly square it k times (`fast=true`) via the tree bipartition scheme to obtain the full-length influence functional.

See [SciPost Phys. Core 7, 063 (2024)].

# Fields
- `algexpan::ExponentialExpansionAlgorithm`: exponential (Prony) expansion algorithm for the bath correlation function.
- `algevo::TimeEvoMPOAlgorithm`: time-evolution algorithm for the differential influence functional (`FirstOrderStepper` or `ComplexStepper`).
- `algmult::DMRGAlgorithm`: MPO multiplication (compression) algorithm.
- `k::Int`: number of tree bipartition iterations, corresponding to 2^k time steps.
- `fast::Bool`: if `true`, use the fast tree bipartition scheme (k multiplications); otherwise use the sequential scheme (2^k-1 multiplications).
- `verbosity::Int`: verbosity level of the output.
"""
struct TranslationInvariantIF{T<:ExponentialExpansionAlgorithm, E<:TimeEvoMPOAlgorithm, M<:DMRGAlgorithm} <: InfluenceFunctionalAlgorithm 
	algexpan::T
	algevo::E
	algmult::M
	# trunc::TruncationDimCutoff
	k::Int
	fast::Bool
	verbosity::Int
end
"""
	TranslationInvariantIF(; algexpan, algevo, algmult, k, fast, verbosity)

Keyword constructor for `TranslationInvariantIF`; all parameters have default values and usually need not be passed explicitly.
"""
TranslationInvariantIF(; algexpan::ExponentialExpansionAlgorithm=OverDeterminedProny(n=15, tol=1.0e-4, verbosity=0), 
						 algevo::TimeEvoMPOAlgorithm=WII(), 
						 algmult::DMRGAlgorithm=DefaultMultAlg,
						 k::Int=5, 
						 fast::Bool=true,
						 verbosity::Int=0) = TranslationInvariantIF(algexpan, algevo, algmult, k, fast, verbosity)

function Base.getproperty(x::TranslationInvariantIF, s::Symbol)
	if s == :trunc
		return x.algmult.trunc
	else
		getfield(x, s)
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
