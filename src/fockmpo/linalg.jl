# function Base.:*(h::FockMPO, psi::FockMPS)
#     @assert !isempty(h)
#     (length(h) == length(psi)) || throw(DimensionMismatch("mpo mps size mismatch"))
#     r = [@tensor tmp[-1 -2; -3 -4 -5] := a[-1, -3, -4, 1] * b[-2, 1, -5] for (a, b) in zip(h.data, psi.data)]
#     return FockMPS([tie(item,(2,1,2)) for item in r], scaling=scaling(psi)*scaling(h))
# end


function Base.:*(hA::FockMPO, hB::FockMPO)
    @assert !isempty(hA)
    (length(hA) == length(hB)) || throw(DimensionMismatch())
    r = [@tensor tmp[1,5,2,3,6,7] := aj[1,2,3,4] * bj[5,4,6,7] for (aj, bj) in zip(hA.data, hB.data)]
    return FockMPO([tie(item, (2,1,2,1)) for item in r])
end

# function integrate(x::FockMPO)
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