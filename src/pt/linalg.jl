

function Base.:*(x::ProcessTensor, y::ProcessTensor)
    @assert !isempty(x)
    (length(x) == length(y)) || throw(DimensionMismatch())
    r = [@tensor tmp[1,5,2,3,6,7] := aj[1,2,3,4] * bj[5,4,6,7] for (aj, bj) in zip(x.data, y.data)]
    return ProcessTensor([tie(item, (2,1,2,1)) for item in r], scaling=scaling(x)*scaling(y))
end



"""
    addition of two MPOs
"""
function Base.:+(hA::ProcessTensor, hB::ProcessTensor)
    @assert !isempty(hA)
    (length(hA) == length(hB)) || throw(DimensionMismatch())
    T = promote_type(scalartype(hA), scalartype(hB))
    L = length(hA)
    scaling_x = scaling(hA)
    scaling_y = scaling(hB)
    (L == 1) && return ProcessTensor([scaling_x * hA[1] + scaling_y * hB[1]])

    r = Vector{Array{T, 4}}(undef, L)
    r[1] = cat(scaling_x * hA[1], scaling_y * hB[1], dims=3)
    r[L] = cat(scaling_x * hA[L], scaling_y * hB[L], dims=1)
    for i in 2:L-1
        r[i] = cat(scaling_x * hA[i], scaling_y * hB[i], dims=(1,3))
    end
    return ProcessTensor(r)
end
# adding mpo with adjoint mpo will return an normal mpo
Base.:-(hA::ProcessTensor, hB::ProcessTensor) = hA + (-1) * hB
Base.:-(h::ProcessTensor) = -1 * h

