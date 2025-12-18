# function integrate(x::ProcessTensor)
#     L = length(x)
#     sca = scaling(x)
#     v = dropdims(x[L], dims=3) * sca
#     for i in L-1:-1:1
#         @tensor tmp[1,5,4] := sca * x[i][1,2,3,4] * v[3,5,2]
#         v = tmp
#     end
#     v = dropdims(v, dims=1)
#     return tr(v)
# end


# function ADT(x::ProcessTensor; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
#     L = length(x)
#     (L > 1) || throw(DimensionMismatch("ProcessTensor size should be larger than 1"))
#     T = scalartype(x)
#     data = Vector{Array{T, 3}}(undef, L+1)
#     d = size(x[1], 2)
#     c = copy_tensor3(T, d)
#     workspace = T[]
#     q, r = tqr!(x[1], (1,4), (2, 3), workspace)
#     data[1] = q
#     for i in 2:L
#         @tensor tmp2[1,7,4,5] := r[1,2,3] * x[i][3,4,5,6] * c[7,2,6]
#         q, r = tqr!(tmp2, (1,2), (3,4), workspace)
#         data[i] = q
#     end
#     data[L+1] = r

#     y = ADT(data)
#     _rightorth!(y, SVD(), trunc, false, verbosity)
#     setscaling!(y, scaling(x)^(L/L+1) * scaling(y))
#     return y
# end


# function copy_tensor3(::Type{T}, d::Int) where {T<:Number}
#     m = zeros(T, d, d, d)
#     for i in 1:d
#         m[i,i,i] = 1
#     end
#     return m
# end