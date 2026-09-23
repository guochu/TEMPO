# 正交化 / 规范化（后端：FiniteMPSAlgorithms，作用在 ProcessTensor 内层的
# CanonicalMPO payload 上，见 fmabackend.jl 与 pt/def.jl）

"""
	leftorth!(h::ProcessTensor; alg::Orthogonalize=Orthogonalize())

Orthogonalize the MPO into left-canonical form, modifying `h` in place and returning it.

# Arguments
- `h::ProcessTensor`: MPO to orthogonalize
- `alg::Orthogonalize`: orthogonalization algorithm configuration; `SVD` without truncation by default

# Returns
`h` itself.
"""
leftorth!(h::ProcessTensor; alg::Orthogonalize = Orthogonalize()) = _leftorth!(h, alg.orth, alg.trunc, alg.normalize, alg.verbosity)
function _leftorth!(h::ProcessTensor, alg::QR, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	v = h.parent
	FiniteMPSAlgorithms._leftorth!(v, alg, trunc, normalize, verbosity)
	return h
end

function _leftorth!(h::ProcessTensor, alg::SVD, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	v = h.parent
	FiniteMPSAlgorithms._leftorth!(v, alg, trunc, normalize, verbosity)
	return h
end

"""
	rightorth!(h::ProcessTensor; alg::Orthogonalize=Orthogonalize())

Orthogonalize the MPO into right-canonical form, modifying `h` in place and returning it.

# Arguments
- `h::ProcessTensor`: MPO to orthogonalize
- `alg::Orthogonalize`: orthogonalization algorithm configuration; `SVD` without truncation by default

# Returns
`h` itself.
"""
rightorth!(h::ProcessTensor; alg::Orthogonalize = Orthogonalize()) = _rightorth!(h, alg.orth, alg.trunc, alg.normalize, alg.verbosity)
function _rightorth!(h::ProcessTensor, alg::QR, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	v = h.parent
	FiniteMPSAlgorithms._rightorth!(v, alg, trunc, normalize, verbosity)
	return h
end

function _rightorth!(h::ProcessTensor, alg::SVD, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	v = h.parent
	FiniteMPSAlgorithms._rightorth!(v, alg, trunc, normalize, verbosity)
	return h
end

canonicalize(psi::ProcessTensor; kwargs...) = canonicalize!(deepcopy(psi); kwargs...)
"""
	canonicalize!(psi::ProcessTensor; alg::Orthogonalize=Orthogonalize(trunc=DefaultITruncation, normalize=false))

Transform the MPO into canonical form, modifying `psi` in place and returning it.

Internally performs a left orthogonalization with `QR` (without truncation) first, then orthogonalizes from right to left using the algorithm specified by `alg` with truncation. Note: enabling normalization (`normalize=true`) is not recommended for `ProcessTensor`.

# Returns
`psi` itself.
"""
function canonicalize!(psi::ProcessTensor; alg::Orthogonalize = Orthogonalize(trunc=DefaultITruncation, normalize=false))
	alg.normalize && @warn "canonicalize with renormalization not recommanded for ProcessTensor"
	v = psi.parent
	FiniteMPSAlgorithms._canonicalize!(v; alg)
	return psi
end
