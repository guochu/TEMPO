

"""
    GenericDecayTerm(a, b, f; middle=_eye(size(a, 1)), coeff=1.0)
    GenericDecayTerm(a, b; middle=_eye(size(a, 1)), f, coeff=1.0)

A generic decaying long-range interaction term of the form `coeff * [â ⊗ f(n) * m̂^⊗n ⊗ b̂]`,
where `f` is a function decaying with the distance `n` (or a pre-sampled vector).
"""
struct GenericDecayTerm{M1<:AbstractMatrix, M<:AbstractMatrix, M2, F, T <: Number} <: AbstractLongRangeTerm
    a::M1
    m::M
    b::M2
    f::F
    coeff::T
end

"""
    GenericDecayTerm(a, b, f; middle=_eye(size(a, 1)), coeff=1.0)

Convenience constructor for `GenericDecayTerm` in which `f` is passed as a positional argument.
"""
GenericDecayTerm(a::AbstractMatrix, b::AbstractMatrix, f; middle::AbstractMatrix = _eye(size(a, 1)), coeff::Number=1.) = GenericDecayTerm(a, middle, b, f, coeff)

"""
    GenericDecayTerm(a, b; middle=_eye(size(a, 1)), f, coeff=1.0)

Keyword constructor for `GenericDecayTerm` in which `f` is passed as a keyword argument.
"""
function GenericDecayTerm(a::AbstractMatrix, b::AbstractMatrix; 
                            middle::AbstractMatrix = _eye(size(a, 1)), f, coeff::Number=1.) 
    GenericDecayTerm(a, middle, b, f, coeff)
end 
TO.scalartype(::Type{GenericDecayTerm{M1, M, M2, F, T}}) where {M1, M, M2, F<:AbstractVector, T} = promote_type(scalartype(M1), scalartype(M), scalartype(M2), eltype(F), T)
TO.scalartype(x::GenericDecayTerm{M1, M, M2, F, T}) where {M1, M, M2, F<:AbstractVector, T} = scalartype(typeof(x))
TO.scalartype(x::GenericDecayTerm{M1, M, M2, F, T}) where {M1, M, M2, F, T} = promote_type(scalartype(M1), scalartype(M), scalartype(M2), T, typeof(x.f(0.)))
Base.adjoint(x::GenericDecayTerm) = GenericDecayTerm(_op_adjoint(x.a, x.m, x.b)..., _conj(x.f), conj(coeff(x)))

_conj(f) = x->conj(f(x))
_conj(f::AbstractVector) = conj(f)

"""
    PowerlawDecayTerm(a::AbstractMatrix, b::AbstractMatrix; α::Number=1., kwargs...)

A power-law decaying long-range interaction term with `f(n) = n^α`. In principle `α` should be negative (otherwise it diverges with distance).
"""
PowerlawDecayTerm(a::AbstractMatrix, b::AbstractMatrix; α::Number=1., kwargs...) = GenericDecayTerm(a, b; f=x->x^α, kwargs...)



# L is the number of sites

"""
    expand_decayterm(x::GenericDecayTerm; len, alg=OverDeterminedProny())

Convert a `GenericDecayTerm` into a list of `ExponentialDecayTerm`s.
When `x.f` is a vector, `len` is ignored; when `x.f` is a function, `len` is the sampling length.
"""
function expand_decayterm(x::GenericDecayTerm{M1, M, M2, F, T}; len::Union{Int, Nothing}=nothing, alg::ExponentialExpansionAlgorithm=OverDeterminedProny()) where {M1, M, M2, F, T}
    if F <: AbstractVector
        xs, lambdas = exponential_expansion(x.f, alg=alg)
        isa(len, Int) && println("key len ignored")
    else
        isa(len, Int) || throw(ArgumentError("key len should be Int when F is not a vector"))
        xs, lambdas = exponential_expansion(x.f, len-1, alg=alg)
    end
    r = ExponentialDecayTerm{M1, M, M2, eltype(lambdas)}[]
    for (c, alpha) in zip(xs, lambdas)
        push!(r, ExponentialDecayTerm(x.a, x.b; middle=x.m, α=alpha, coeff=c * coeff(x)))
    end
    return r
end
