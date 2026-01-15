# orthogonalize mps to be left-canonical or right-canonical

leftorth!(h::ProcessTensor; alg::Orthogonalize = Orthogonalize()) = _leftorth!(h, alg.orth, alg.trunc, alg.normalize, alg.verbosity)
function _leftorth!(psi::ProcessTensor, alg::QR, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	L = length(psi)
	workspace = scalartype(psi)[]
	for i in 1:L-1
		q, r = tqr!(psi[i], (1, 2, 4), (3,), workspace)
		_renormalize!(psi, r, normalize)
		psi[i] = permute(q, (1,2,4,3))
		psi[i+1] = @tensor tmp[1,3,4,5] := r[1,2] * psi[i+1][2,3,4,5]
	end
	_renormalize!(psi, psi[L], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

function _leftorth!(psi::ProcessTensor, alg::SVD, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	L = length(psi)
	workspace = scalartype(psi)[]
	maxerr = 0.
	for i in 1:L-1
		u, s, v, err = tsvd!(psi[i], (1,2,4), (3), workspace, trunc=trunc)
		nr = _renormalize!(psi, s, normalize)
		rerror = sqrt(err * err / (nr * nr + err * err))
		(verbosity > 1) && println("SVD truncerror at bond $(i): ", rerror)
		psi[i] = permute(u, (1,2,4,3))
		v = Diagonal(s) * v
		psi[i+1] = @tensor tmp[1,3,4,5] := v[1,2] * psi[i+1][2,3,4,5]
		psi.s[i+1] = s
		maxerr = max(maxerr, rerror)
	end
	(verbosity > 0) && println("Max SVD truncerror in leftorth: ", maxerr)
	_renormalize!(psi, psi[L], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

rightorth!(h::ProcessTensor; alg::Orthogonalize = Orthogonalize()) = _rightorth!(h, alg.orth, alg.trunc, alg.normalize, alg.verbosity)
function _rightorth!(psi::ProcessTensor, alg::QR, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	L = length(psi)
	workspace = scalartype(psi)[]
	for i in L:-1:2
		l, q = tlq!(psi[i], (1,), (2, 3, 4), workspace)
		_renormalize!(psi, l, normalize)
		psi[i] = q
		psi[i-1] = @tensor tmp[1,2,5,4] := psi[i-1][1,2,3,4] * l[3,5]
	end
	_renormalize!(psi, psi[1], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end
function _rightorth!(psi::ProcessTensor, alg::SVD, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	L = length(psi)
	workspace = scalartype(psi)[]
	maxerr = 0.
	for i in L:-1:2
		u, s, v, err = tsvd!(psi[i], (1,), (2, 3, 4), workspace, trunc=trunc)
		psi[i] = v
		nr = _renormalize!(psi, s, normalize)
		rerror = sqrt(err * err / (nr * nr + err * err))
		(verbosity > 1) && println("SVD truncerror at bond $(i): ", rerror)
		u = u * Diagonal(s)
		psi[i-1] = @tensor tmp[1,2,5,4] := psi[i-1][1,2,3,4] * u[3,5]
		psi.s[i] = s
		maxerr = max(maxerr, err)
	end
	(verbosity > 0) && println("Max SVD truncerror in rightorth: ", maxerr)
	_renormalize!(psi, psi[1], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

canonicalize(psi::ProcessTensor; kwargs...) = canonicalize!(deepcopy(psi); kwargs...)
function canonicalize!(psi::ProcessTensor; alg::Orthogonalize = Orthogonalize(trunc=DefaultTruncation, normalize=false))
	alg.normalize && @warn "canonicalize with renormalization not recommanded for ProcessTensor"
	L = length(psi)
	_leftorth!(psi, QR(), NoTruncation(), alg.normalize, alg.verbosity)
	_rightorth!(psi, alg.orth, alg.trunc, alg.normalize, alg.verbosity)
	return psi
end