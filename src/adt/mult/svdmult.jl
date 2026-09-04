"""
	mult!(x::ADT, y::ADT; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)

Compute the product (tensor contraction) of two MPS in place on `x`, compress the result, and return `x`.

Algorithm: accumulate left-to-right QR decompositions, then orthogonalize right-to-left with `SVD` and truncate according to `trunc`; the resulting scaling factor is `scaling(x) * scaling(y)`.

# Arguments
- `x::ADT`: first MPS, also used as the output storage
- `y::ADT`: second MPS, must have the same length as `x`
- `trunc`: SVD truncation scheme (e.g., `truncdim(D)`, `trunccutoff(ε)`)
- `verbosity::Int`: verbosity level

# Returns
`x` itself (its content is overwritten by the product).
"""
function mult!(x::ADT, y::ADT; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
    (length(x) == length(y)) || throw(DimensionMismatch())
    T = promote_type(scalartype(x), scalartype(y))
    L = length(x)
    z = Vector{Array{T, 3}}(undef, L)
    tmp = tie(n_fuse(_mult_site_n(x[1], y[1]), 3), (2,1,1,1))
    q, r = leftorth!(tmp, (1, 2), (3,4))
    x[1] = q
    for i in 2:L-1
        @tensor tmp[1,4,6,5,7] := r[1,2,3] * x[i][2,4,5] * y[i][3,6,7]
        tmp2 = n_fuse(tmp, 2)
        q, r = leftorth!(tmp2, (1, 2), (3,4))
        x[i] = q
        _renormalize!(x, r, false)
    end
    @tensor tmp[1,4,6,5,7] := r[1,2,3] * x[L][2,4,5] * y[L][3,6,7]
    x[end] = tie(n_fuse(tmp, 2), (1,2,1))
    _rightorth!(x, SVD(), trunc, false, verbosity)
    setscaling!(x, scaling(x) * scaling(y))
    return x
end
"""
	mult(x::ADT, y::ADT; kwargs...)

Compute the product of two MPS and compress it, returning a new `ADT` (inputs are not modified).

Equivalent to `mult!(copy(x), y; kwargs...)`; arguments as in `mult!`.

# Returns
The product `ADT`.
"""
mult(x::ADT, y::ADT; kwargs...) = mult!(copy(x), y; kwargs...)