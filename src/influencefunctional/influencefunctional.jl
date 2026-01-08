abstract type InfluenceFunctionalAlgorithm end
"""
	struct PartialIF

Build the IF as the product of partial MPOs, each with D=2
see [SciPost Phys. Core 7, 063 (2024)]
"""
struct PartialIF <: InfluenceFunctionalAlgorithm 
	trunc::TruncationDimCutoff
end
PartialIF(; trunc::TruncationDimCutoff=DefaultITruncation) = PartialIF(trunc)

"""
	struct TranslationInvariantIF

Build the IF as a translational variant MPO
see [SciPost Phys. Core 7, 063 (2024)]
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
TranslationInvariantIF(; algexpan::ExponentialExpansionAlgorithm=PronyExpansion(n=15, tol=1.0e-4, verbosity=0), 
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




abstract type HybridizationStyle end


struct AdditiveHyb <: HybridizationStyle
	op::Vector{Float64}
end

AdditiveHyb(x::AbstractVector{<:Real}) = AdditiveHyb(float(x))

phydim(b::AdditiveHyb) = length(b.op)

TO.scalartype(::Type{AdditiveHyb}) = Float64

abstract type GeneralHybStyle <: HybridizationStyle end

struct NonAdditiveHyb{T<:Number} <: GeneralHybStyle
	op::Matrix{T}

function NonAdditiveHyb{T}(op::AbstractMatrix) where {T<:Number}
	(size(op, 1) == size(op, 2)) || throw(ArgumentError("square matrix expected"))
	new{T}(convert(Matrix{T}, op))
end
end
NonAdditiveHyb(a::AbstractMatrix{T}) where {T<:Number} = NonAdditiveHyb{T}(a)

phydim(b::NonAdditiveHyb) = size(b.op, 1)

TO.scalartype(::Type{NonAdditiveHyb{T}}) where T = T
pairop(b::NonAdditiveHyb) = b.op, b.op

"""
	struct NonDiagonalHyb{T<:Number}

The Impurity-bath coupling is op*a + op'*a'
"""
struct NonDiagonalHyb{T<:Number} <: GeneralHybStyle
	op::Matrix{T}

function NonDiagonalHyb{T}(op::AbstractMatrix) where {T<:Number}
	(size(op, 1) == size(op, 2)) || throw(ArgumentError("square matrix expected"))
	new{T}(convert(Matrix{T}, op))
end
end
NonDiagonalHyb(a::AbstractMatrix{T}) where {T<:Number} = NonDiagonalHyb{T}(a)

phydim(b::NonDiagonalHyb) = size(b.op, 1)

TO.scalartype(::Type{NonDiagonalHyb{T}}) where T = T
pairop(b::NonDiagonalHyb) = b.op, b.op'

include("partialif/partialif.jl")
include("ptpartialif/ptpartialif.jl")
include("ttiif/ttiif.jl")
