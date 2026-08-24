"""
	mult!(x::ProcessTensor, y::ProcessTensor; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)

Compute the product (tensor contraction) of two MPO in place on `x`, compress the result, and return `x`.

Algorithm: accumulate left-to-right QR decompositions, then orthogonalize right-to-left with `SVD` and truncate according to `trunc`; the resulting scaling factor is `scaling(x) * scaling(y)`.

# Arguments
- `x::ProcessTensor`: first MPO, also used as the output storage
- `y::ProcessTensor`: second MPO, must have the same length as `x`
- `trunc`: SVD truncation scheme (e.g., `truncdim(D)`, `trunccutoff(ε)`)
- `verbosity::Int`: verbosity level

# Returns
`x` itself (its content is overwritten by the product).
"""
function mult!(x::ProcessTensor, y::ProcessTensor; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
    (length(x) == length(y)) || throw(DimensionMismatch())
    T = promote_type(scalartype(x), scalartype(y))
    L = length(x)
    z = Vector{Array{T, 3}}(undef, L)
    tmp = tie(_mult_mpo_sitetensor(x[1], y[1]), (2,1,1,1,1))
    workspace = scalartype(tmp)[]
    q, r = tqr!(tmp, (1, 2, 5), (3,4), workspace)
    x[1] = permute(q, (1,2,4,3))
    for i in 2:L-1
        @tensor tmp[1,4,5,7,8] := r[1,2,3] * x[i][2,4,5,6] * y[i][3,6,7,8]
        q, r = tqr!(tmp, (1, 2, 5), (3,4), workspace)
        x[i] = permute(q, (1,2,4,3))
        _renormalize!(x, r, false)
    end
    @tensor tmp[1,4,5,7,8] := r[1,2,3] * x[L][2,4,5,6] * y[L][3,6,7,8]
    x[end] = tie(tmp, (1,1,2,1))
    _rightorth!(x, SVD(), trunc, false, verbosity)
    setscaling!(x, scaling(x) * scaling(y))
    return x
end
"""
	mult(x::ProcessTensor, y::ProcessTensor; kwargs...)

Compute the product of two MPO and compress it, returning a new `ProcessTensor` (inputs are not modified).

Equivalent to `mult!(copy(x), y; kwargs...)`; arguments as in `mult!`.

# Returns
The product `ProcessTensor`.
"""
mult(x::ProcessTensor, y::ProcessTensor; kwargs...) = mult!(copy(x), y; kwargs...)


function _mult_mpo_sitetensor(xj::DenseMPOTensor, yj::DenseMPOTensor)
    @tensor tmp[1,5,2,3,6,7] := xj[1,2,3,4] * yj[5,4,6,7]
    return tmp
end