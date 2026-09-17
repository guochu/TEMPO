

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


function easy_swap!(x::ProcessTensor, bond::Int; trunc::TruncationScheme=DefaultTruncation)
	x[bond], x.s[bond+1], x[bond+1] = _swap_gate(x.s[bond], x[bond], x.s[bond+1], x[bond+1], trunc=trunc)
	return x
end

# Hastings-style swap gate (following GTEMPO): the bond Schmidt values `svectorj1`
# stored on the left of the swapped pair are contracted into the two-site block,
# which is then re-decomposed; the renewed bond spectrum and the right-canonical
# factor are written back to `svectorj2` and the second site tensor.
# Site tensor layout: (aL, pout, aR, pin).
function _swap_gate(svectorj1::Vector, m1::DenseMPOTensor, svectorj2::Vector, m2::DenseMPOTensor; trunc::TruncationScheme)
	sv1 = Diagonal(svectorj1)
	local twositemps
	@tensor twositemps[a, b, c, d, e, f] := m1[a, b, 2, c] * m2[2, d, f, e]
	local twositemps1
	@tensor twositemps1[a, b, c, d, e, f] := sv1[a, 1] * twositemps[1, b, c, d, e, f]
	u, s, v = tsvd!(twositemps1, (1, 2, 3), (4, 5, 6); trunc=trunc)
	local u2
	@tensor u2[a, b, c, d] := twositemps[a, b, c, 1, 2, 3] * conj(v[d, 1, 2, 3])
	return permute(u2, (1, 2, 4, 3)), s, permute(v, (1, 2, 4, 3))
end


function _permute!(x::ProcessTensor, perm::Vector{Int}; trunc::TruncationScheme=DefaultIntegrationTruncation)
	@assert length(x) == length(perm)
	if svectors_uninitialized(x)
		canonicalize!(x, alg=Orthogonalize(trunc=trunc, normalize=false))
	end
	p = permutation2swaps(perm)
	for i in p
		easy_swap!(x, i, trunc=trunc)
	end
	return x
end
permute!(x::ProcessTensor, perm::Vector; kwargs...) = _permute!(x, perm; kwargs...)
permute(x::ProcessTensor, perm::Vector{Int}; kwargs...) = permute!(deepcopy(x), perm; kwargs...)

