# TEMPO 特有的长程衰减项（long-range decay terms）
#
# GenericDecayTerm / PowerlawDecayTerm 描述任意衰减函数的长程相互作用，
# `expand_decayterm` 经指数展开（Prony）打包为 FiniteMPSAlgorithms 的
# `ExpDecayOpSum`，再由 `SchurMPOTensor(::ExpDecayOpSum, hloc)` 构造
# Schur 形式的 MPO 站点张量（每个展开参数一个内部通道）。

# ---- 衰减项的公共接口 ----

abstract type AbstractLongRangeTerm end

space_l(x::AbstractLongRangeTerm) = isa(x.a, AbstractMatrix) ? size(x.a, 1) : 1
space_r(x::AbstractLongRangeTerm) = isa(x.b, AbstractMatrix) ? size(x.b, 2) : 1
coeff(x::AbstractLongRangeTerm) = x.coeff

# ---- GenericDecayTerm / PowerlawDecayTerm ----

"""
    GenericDecayTerm(a, b, f; middle=isometry(size(a, 1)), coeff=1.0)
    GenericDecayTerm(a, b; middle=isometry(size(a, 1)), f, coeff=1.0)

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
    GenericDecayTerm(a, b, f; middle=isometry(size(a, 1)), coeff=1.0)

Convenience constructor for `GenericDecayTerm` in which `f` is passed as a positional argument.
"""
GenericDecayTerm(a::AbstractMatrix, b::AbstractMatrix, f; middle::AbstractMatrix = isometry(size(a, 1)), coeff::Number=1.) = GenericDecayTerm(a, middle, b, f, coeff)

"""
    GenericDecayTerm(a, b; middle=isometry(size(a, 1)), f, coeff=1.0)

Keyword constructor for `GenericDecayTerm` in which `f` is passed as a keyword argument.
"""
function GenericDecayTerm(a::AbstractMatrix, b::AbstractMatrix;
                            middle::AbstractMatrix = isometry(size(a, 1)), f, coeff::Number=1.)
    GenericDecayTerm(a, middle, b, f, coeff)
end
TO.scalartype(::Type{GenericDecayTerm{M1, M, M2, F, T}}) where {M1, M, M2, F<:AbstractVector, T} = promote_type(scalartype(M1), scalartype(M), scalartype(M2), eltype(F), T)
TO.scalartype(x::GenericDecayTerm{M1, M, M2, F, T}) where {M1, M, M2, F<:AbstractVector, T} = scalartype(typeof(x))
TO.scalartype(x::GenericDecayTerm{M1, M, M2, F, T}) where {M1, M, M2, F, T} = promote_type(scalartype(M1), scalartype(M), scalartype(M2), T, typeof(x.f(0.)))
Base.adjoint(x::GenericDecayTerm) = GenericDecayTerm(_op_adjoint(x.a, x.m, x.b)..., _conj(x.f), conj(coeff(x)))

_op_adjoint(a::AbstractMatrix, m::AbstractMatrix, b::AbstractMatrix) = (a', m', b')
_conj(f) = x->conj(f(x))
_conj(f::AbstractVector) = conj(f)

"""
    PowerlawDecayTerm(a::AbstractMatrix, b::AbstractMatrix; α::Number=1., kwargs...)

A power-law decaying long-range interaction term with `f(n) = n^α`. In principle `α` should be negative (otherwise it diverges with distance).
"""
PowerlawDecayTerm(a::AbstractMatrix, b::AbstractMatrix; α::Number=1., kwargs...) = GenericDecayTerm(a, b; f=x->x^α, kwargs...)

# ---- 指数展开（Prony）----

# L is the number of sites

"""
    expand_decayterm(x::GenericDecayTerm; len, alg=OverDeterminedProny())

Convert a `GenericDecayTerm` into an `ExpDecayOpSum`: the exponential expansion
coefficients (times the term's `coeff`) become the strengths `αs`, the exponents the
decay factors `λs`. The sum converts directly to a `SchurMPOTensor` (one internal
channel per expansion parameter). When `x.f` is a vector, `len` is ignored; when
`x.f` is a function, `len` is the sampling length.
"""
function expand_decayterm(x::GenericDecayTerm{M1, M, M2, F, T}; len::Union{Int, Nothing}=nothing, alg::ExponentialExpansionAlgorithm=OverDeterminedProny()) where {M1, M, M2, F, T}
    if F <: AbstractVector
        xs, lambdas = exponential_expansion(x.f, alg=alg)
        isa(len, Int) && println("key len ignored")
    else
        isa(len, Int) || throw(ArgumentError("key len should be Int when F is not a vector"))
        xs, lambdas = exponential_expansion(x.f, len-1, alg=alg)
    end
    return ExpDecayOpSum(x.a, x.m, x.b, [c * coeff(x) for c in xs], lambdas)
end
