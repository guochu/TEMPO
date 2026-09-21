

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


function swap!(x::ProcessTensor, bond::Int; trunc::TruncationScheme=DefaultITruncation)
	x[bond], x.s[bond+1], x[bond+1] = _swap_gate(x[bond], x[bond+1], trunc=trunc)
	return x
end

# Swap gate: builds the two-site block `x[bond] · x[bond+1]` (the PT contraction
# convention carries the bond spectra inside the site tensors, so no explicit
# spectrum is contracted here) with the conjugate pairs already swapped, i.e. in
# the order (aL, pout2, pin2, pout1, pin1, aR), and re-decomposes it such that
# the left factor carries the conjugate pair of site bond+1 and the right factor
# the one of site bond. The renewed bond spectrum is absorbed into the left
# factor (its norm equals the bond spectrum recorded in `x.s[bond+1]`), so the
# chain contraction is preserved while the sites are exchanged.
function _swap_gate(m1::DenseMPOTensor, m2::DenseMPOTensor; trunc::TruncationScheme)
	@tensor block[a, e, f, b, c, d] := m1[a, b, k, c] * m2[k, e, d, f]
	u, s, v = tsvd!(block, (1, 2, 3), (4, 5, 6); trunc=trunc)
	u = u .* reshape(Vector(s), 1, 1, 1, :)
	return permute(u, (1, 2, 4, 3)), s, permute(v, (1, 2, 4, 3))
end


function _permute!(x::ProcessTensor, perm::Vector{Int}; trunc::TruncationScheme=DefaultKTruncation)
	@assert length(x) == length(perm)
	if svectors_uninitialized(x)
		canonicalize!(x, alg=Orthogonalize(trunc=trunc, normalize=false))
	end
	p = permutation2swaps(perm)
	for i in p
		swap!(x, i, trunc=trunc)
	end
	return x
end
permute!(x::ProcessTensor, perm::Vector; kwargs...) = _permute!(x, perm; kwargs...)
permute(x::ProcessTensor, perm::Vector{Int}; kwargs...) = permute!(deepcopy(x), perm; kwargs...)

