

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
	local twositemps
	# site tensor layout: (aL, pout, aR, pin)
	@tensor twositemps[a, b, c, d, e, f] := x[bond][a, b, 2, c] * x[bond+1][2, d, f, e]
	u, s, v = tsvd!(twositemps, (1, 2, 3), (4, 5, 6); trunc=trunc)
	x[bond] = permute(u .* reshape(s, 1, 1, 1, :), (1, 2, 4, 3))
	x[bond+1] = permute(v, (1, 2, 4, 3))
	x.s[bond+1] = s
	return x
end

function naive_swap!(x::ProcessTensor, bond::Int; trunc::TruncationScheme=DefaultTruncation)
	x[bond], x[bond+1] = _swap_gate(x[bond], x[bond+1], trunc=trunc)
	return x
end

# swap gate of two adjacent site tensors, absorbing the Schmidt spectrum into the left factor
function _swap_gate(m1::DenseMPOTensor, m2::DenseMPOTensor; trunc::TruncationScheme)
	@tensor twositemps[1, 2, 4, 5, 3, 6] := m1[1, 2, 3, 4] * m2[4, 5, 6, 7]
	u, s, v = tsvd!(twositemps, (1, 2, 5), (3, 4, 6); trunc=trunc)
	return u .* reshape(s, 1, 1, 1, :), permute(v, (1, 2, 4, 3))
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
	# the swaps reshuffle the entanglement across bonds: re-establish the
	# mixed-canonical form (and the bond Schmidt records) at the end
	canonicalize!(x, alg=Orthogonalize(trunc=trunc, normalize=false))
	return x
end
permute!(x::ProcessTensor, perm::Vector; kwargs...) = _permute!(x, perm; kwargs...)
permute(x::ProcessTensor, perm::Vector{Int}; kwargs...) = permute!(deepcopy(x), perm; kwargs...)

function naive_permute!(x::ProcessTensor, perm::Vector{Int}; trunc::TruncationScheme=DefaultIntegrationTruncation)
	@assert length(x) == length(perm)
	p = permutation2swaps(perm)
	for i in p
		naive_swap!(x, i, trunc=trunc)
	end
	return x
end
naive_permute(x::ProcessTensor, perm::Vector{Int}; kwargs...) = naive_permute!(copy(x), perm; kwargs...)

