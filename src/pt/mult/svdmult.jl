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
mult(x::ProcessTensor, y::ProcessTensor; kwargs...) = mult!(copy(x), y; kwargs...)


function _mult_mpo_sitetensor(xj::DenseMPOTensor, yj::DenseMPOTensor)
    @tensor tmp[1,5,2,3,6,7] := xj[1,2,3,4] * yj[5,4,6,7]
    return tmp
end