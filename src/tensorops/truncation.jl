"""
    TruncationScheme

Abstract type for tensor truncation schemes, specifying how singular values are truncated (by dimension, by truncation error, or both).
"""
abstract type TruncationScheme end


"""
    NoTruncation()

A truncation scheme that performs no truncation, keeping all singular values.
"""
struct NoTruncation <: TruncationScheme end

struct TruncateDim <: TruncationScheme
	D::Int
end
TruncateDim(;D::Int) = TruncateDim(D)
"""
    truncdim(d::Int)

Construct a dimension-truncated scheme `TruncateDim(d)` that keeps only the `d` largest singular values.
"""
truncdim(d::Int) = TruncateDim(d)
"""
    truncdim(; D::Int)

Keyword form of the dimension-truncation scheme constructor, equivalent to `truncdim(D)`.
"""
truncdim(; D::Int) = truncdim(D)

struct TruncateCutoff <: TruncationScheme
	ϵ::Float64
end
TruncateCutoff(;ϵ::Real) = TruncateCutoff(convert(Float64, ϵ))
"""
    trunccutoff(; ϵ::Real)

Construct a cutoff-truncated scheme `TruncateCutoff(ϵ)` that discards singular values with relative norm below `ϵ`.
"""
trunccutoff(; ϵ::Real) = TruncateCutoff(ϵ)

# reserve at least add_back singular values
"""
    TruncationDimCutoff(D, ϵ, add_back=0)
    TruncationDimCutoff(; D, ϵ, add_back=0)

A truncation scheme combining dimension and cutoff: the truncation point is first determined by the relative norm `ϵ`,
then capped at `D` singular values while keeping at least `add_back` of them.
"""
struct TruncationDimCutoff <: TruncationScheme
    D::Int
    ϵ::Float64
    add_back::Int
end
TruncationDimCutoff(;D::Int, ϵ::Real, add_back::Int=0) = TruncationDimCutoff(D, float(ϵ), min(add_back, D))
"""
    truncdimcutoff(D, ϵ, add_back=0)

Positional-argument convenience constructor for `TruncationDimCutoff`, equivalent to `TruncationDimCutoff(D, ϵ, add_back)`.
"""
truncdimcutoff(D::Int, epsilon::Real; add_back::Int=0) = TruncationDimCutoff(D, epsilon, min(add_back, D))
"""
    truncdimcutoff(; D, ϵ, add_back=0)

Keyword convenience constructor for `TruncationDimCutoff`, equivalent to `TruncationDimCutoff(D, ϵ, add_back)`.
"""
truncdimcutoff(; D::Int, ϵ::Real, add_back::Int=0) = TruncationDimCutoff(D, float(ϵ), min(add_back, D))


_truncate!(v::AbstractVector{<:Real}, trunc::NoTruncation, p::Real=2) = v, 0.

function _truncate!(v::AbstractVector{<:Real}, trunc::TruncateDim, p::Real=2)
	dtrunc = min(length(v), trunc.D)
	truncerr = norm(view(v, dtrunc+1:length(v)), p)
	resize!(v, dtrunc)
	return v, truncerr
end

function _truncate!(v::AbstractVector{<:Real}, trunc::TruncateCutoff, p::Real=2)
	sca = norm(v, p)
	dtrunc = findlast(Base.Fix2(>, sca * trunc.ϵ), v)
	if isnothing(dtrunc)
		dtrunc = 0
	end
	return _truncate!(v, TruncateDim(dtrunc), p)
end

function _truncate!(v::AbstractVector{<:Real}, trunc::TruncationDimCutoff, p::Real=2)
	sca = norm(v, p)
	dtrunc = findlast(Base.Fix2(>, sca * trunc.ϵ), v)
	dtrunc = isnothing(dtrunc) ? 0 : dtrunc
	dtrunc = max(dtrunc, trunc.add_back)   # keep at least add_back singular values
	dtrunc = min(dtrunc, trunc.D)          # but never more than D
	v, err = _truncate!(v, TruncateDim(dtrunc), p)
	return v, err / sca
end


# const DefaultTruncation = truncdimcutoff(D=200, ϵ=1.0e-6)
