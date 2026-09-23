# 正交化 / 规范化（后端：FiniteMPSAlgorithms）
#
# `Orthogonalize` 来自 FiniteMPSAlgorithms；实际的 QR/SVD sweep 由
# FiniteMPSAlgorithms 在 ADT 内层的 CanonicalMPS payload 上完成
# （见 fmabackend.jl 与 adt/def.jl）。

"""
	leftorth!(psi::ADT; alg::Orthogonalize=Orthogonalize())

Orthogonalize the MPS into left-canonical form, modifying `psi` in place and returning it.

# Arguments
- `psi::ADT`: MPS to orthogonalize
- `alg::Orthogonalize`: orthogonalization algorithm configuration; `SVD` without truncation by default

# Returns
`psi` itself.
"""
leftorth!(psi::ADT; alg::Orthogonalize = Orthogonalize()) = _leftorth!(psi, alg.orth, alg.trunc, alg.normalize, alg.verbosity)
function _leftorth!(psi::ADT, alg::QR, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	v = psi.parent
	FiniteMPSAlgorithms._leftorth!(v, alg, trunc, normalize, verbosity)
	return psi
end

function _leftorth!(psi::ADT, alg::SVD, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	v = psi.parent
	FiniteMPSAlgorithms._leftorth!(v, alg, trunc, normalize, verbosity)
	return psi
end

"""
	rightorth!(psi::ADT; alg::Orthogonalize=Orthogonalize())

Orthogonalize the MPS into right-canonical form, modifying `psi` in place and returning it.

# Arguments
- `psi::ADT`: MPS to orthogonalize
- `alg::Orthogonalize`: orthogonalization algorithm configuration; `SVD` without truncation by default

# Returns
`psi` itself.
"""
rightorth!(psi::ADT; alg::Orthogonalize = Orthogonalize()) = _rightorth!(psi, alg.orth, alg.trunc, alg.normalize, alg.verbosity)
function _rightorth!(psi::ADT, alg::QR, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	v = psi.parent
	FiniteMPSAlgorithms._rightorth!(v, alg, trunc, normalize, verbosity)
	return psi
end

function _rightorth!(psi::ADT, alg::SVD, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	v = psi.parent
	FiniteMPSAlgorithms._rightorth!(v, alg, trunc, normalize, verbosity)
	return psi
end

canonicalize(psi::ADT; kwargs...) = canonicalize!(deepcopy(psi); kwargs...)
"""
	canonicalize!(psi::ADT; alg::Orthogonalize=Orthogonalize(trunc=DefaultITruncation, normalize=false))

Transform the MPS into canonical form, modifying `psi` in place and returning it.

Internally performs a left orthogonalization with `QR` (without truncation) first, then orthogonalizes from right to left using the algorithm specified by `alg` with truncation. Note: enabling normalization (`normalize=true`) is not recommended for `ADT`.

# Arguments
- `psi::ADT`: MPS to orthogonalize
- `alg::Orthogonalize`: algorithm configuration for the right-orthogonalization stage

# Returns
`psi` itself.
"""
function canonicalize!(psi::ADT; alg::Orthogonalize = Orthogonalize(trunc=DefaultITruncation, normalize=false))
	alg.normalize && @warn "canonicalize with renormalization not recommanded for ADT"
	v = psi.parent
	FiniteMPSAlgorithms._canonicalize!(v; alg)
	return psi
end

function _rescaling!(psi, n::Real)
	L = length(psi)
	scale1 = n^(1/L)
	setscaling!(psi, scaling(psi) * scale1)
	return psi
end
function _rescaling!(psi)
	nrm1 = norm(psi[1])
	psi[1] = rmul!(psi[1], 1/nrm1)
	return _rescaling!(psi, nrm1)
end
