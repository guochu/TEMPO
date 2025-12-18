# # # the convention here is different from DMRG!!!
LinearAlgebra.dot(psiA::Dense1DTN, psiB::Dense1DTN) = _dot(psiA, psiB) * (scaling(psiA) * scaling(psiB))^length(psiA)
function LinearAlgebra.norm(psi::Dense1DTN) 
	a = real(_dot(psi, psi))
    a = (abs(a) >= 1.0e-14) ? a : zero(a)
	return sqrt(a) * scaling(psi)^(length(psi))
end
distance(a::Dense1DTN, b::Dense1DTN) = _distance(a, b)
distance2(a::Dense1DTN, b::Dense1DTN) = _distance2(a, b)


function LinearAlgebra.lmul!(f::Number, psi::Dense1DTN)
    if !isempty(psi)
        psi[1] *= f
    end
    _renormalize!(psi, psi[1], false)
    return psi
end

Base.:*(psi::Dense1DTN, f::Number) = lmul!(f, copy(psi))
Base.:*(f::Number, psi::Dense1DTN) = psi * f
Base.:/(psi::Dense1DTN, f::Number) = psi * (1/f)
Base.:(-)(psi::Dense1DTN) = (-1) * psi

function _dot(psiA::Dense1DTN, psiB::Dense1DTN) 
    (length(psiA) == length(psiB)) || throw(ArgumentError("dimension mismatch"))
    hold = l_LL(psiA, psiB)
    for i in 1:length(psiA)
        hold = updateleft(hold, psiA[i], psiB[i])
    end
    return tr(hold)

end


# the reuslt is also a GrassmannMPS
function Base.:*(x::ADT, y::ADT)
    (length(x) == length(y)) || throw(DimensionMismatch())
    r = [n_fuse(_mult_site_n(x[i], y[i]), 3) for i in 1:length(x)]
    return ADT([tie(rj, (2,1,2)) for rj in r], scaling=scaling(x)*scaling(y))
end

function Base.:+(x::ADT, y::ADT) 
    (length(x) == length(y)) || throw(DimensionMismatch())
    @assert !isempty(x)
    scaling_x = scaling(x)
    scaling_y = scaling(y)
    (length(x) == 1) && return ADT([scaling_x * x[1] + scaling_y * y[1]])

    L = length(x)
    T = promote_type(scalartype(x), scalartype(y))
    r = Vector{Array{T, 3}}(undef, L)
    r[1] = cat(scaling_x*x[1], scaling_y*y[1], dims=3)
    r[L] = cat(scaling_x*x[L], scaling_y*y[L], dims=1)
    for i in 2:L-1
        r[i] = cat(scaling_x*x[i], scaling_y*y[i], dims=(1,3))
    end
    return ADT(r)
end
Base.:-(x::ADT, y::ADT) = x + (-y)


function _permute!(x::ADT, perm::Vector{Int}; trunc::TruncationScheme=DefaultIntegrationTruncation)
    @assert length(x) == length(perm)
    if svectors_uninitialized(x)
        canonicalize!(x, alg=Orthogonalize(trunc=trunc, normalize=false))
    end
    p = CoxeterDecomposition(Permutation(perm))
    for i in p.terms
        easy_swap!(x, i, trunc=trunc)
    end
    return x
end
permute!(x::ADT, perm::Vector; kwargs...) = _permute!(x, perm; kwargs...)
permute(x::ADT, perm::Vector{Int}; kwargs...) = permute!(deepcopy(x), perm; kwargs...)

function naive_permute!(x::ADT, perm::Vector{Int}; trunc::TruncationScheme=DefaultIntegrationTruncation)
    @assert length(x) == length(perm)
    p = CoxeterDecomposition(Permutation(perm))
    for i in p.terms
        naive_swap!(x, i, trunc=trunc)
    end
    return x
end
naive_permute(x::ADT, perm::Vector{Int}; kwargs...) = naive_permute!(copy(x), perm; kwargs...)

function _mult_site_n(xj::DenseMPSTensor, yj::DenseMPSTensor)
    @tensor r[1,4,2,5;3,6] := xj[1,2,3] * yj[4,5,6]
    return r
end

