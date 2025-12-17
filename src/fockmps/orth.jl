# orthogonalize mps to be left-canonical or right-canonical
abstract type MatrixProductOrthogonalAlgorithm  end

"""
	struct MatrixProductOrthogonalize{A<:Union{QR, SVD}, T<:TruncationScheme}
"""
struct Orthogonalize{A<:Union{QR, SVD}, T<:TruncationScheme} <: MatrixProductOrthogonalAlgorithm
	orth::A
	trunc::T
	normalize::Bool
	verbosity::Int
end
Orthogonalize(a::Union{QR, SVD}, trunc::TruncationScheme; normalize::Bool=false, verbosity::Int=0) = Orthogonalize(a, trunc, normalize, verbosity)
Orthogonalize(a::Union{QR, SVD}; trunc::TruncationScheme=NoTruncation(), normalize::Bool=false, verbosity::Int=0) = Orthogonalize(a, trunc, normalize, verbosity)
Orthogonalize(; alg::Union{QR, SVD} = SVD(), trunc::TruncationScheme=NoTruncation(), normalize::Bool=false, verbosity::Int=0) = Orthogonalize(alg, trunc, normalize, verbosity)


leftorth!(psi::FockMPS; alg::Orthogonalize = Orthogonalize()) = _leftorth!(psi, alg.orth, alg.trunc, alg.normalize, alg.verbosity)
function _leftorth!(psi::FockMPS, alg::QR, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	L = length(psi)
	for i in 1:L-1
		q, r = tqr!(psi[i], (1, 2), (3,))
		psi[i] = q
		_renormalize!(psi, r, normalize)
		@tensor tmp[1 3; 4] := r[1,2] * psi[i+1][2,3,4]
		psi[i + 1] = tmp
	end
	_renormalize!(psi, psi[L], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

function _leftorth!(psi::FockMPS, alg::SVD, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	L = length(psi)
	# errs = Float64[]
	maxerr = 0.
	for i in 1:L-1
		u, s, v, err = tsvd!(psi[i], (1, 2), (3,), trunc=trunc)
		nr = _renormalize!(psi, s, normalize)
		rerror = sqrt(err * err / (nr * nr + err * err))
		(verbosity > 1) && println("SVD truncerror at bond $(i): ", rerror)
		psi[i] = u
		v2 = Diagonal(s) * v
		@tensor tmp[-1 -2; -3] := v2[-1, 1] * psi[i+1][1,-2,-3]
		psi[i+1] = tmp
		psi.s[i+1] = s
		# push!(errs, err)
		maxerr = max(maxerr, rerror)
	end
	(verbosity > 0) && println("Max SVD truncerror in leftorth: ", maxerr)
	_renormalize!(psi, psi[L], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

rightorth!(psi::FockMPS; alg::Orthogonalize = Orthogonalize()) = _rightorth!(psi, alg.orth, alg.trunc, alg.normalize, alg.verbosity)
function _rightorth!(psi::FockMPS, alg::QR, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	L = length(psi)
	for i in L:-1:2
		l, q = tlq!(psi[i], (1,), (2, 3))
		psi[i] = q
		_renormalize!(psi, l, normalize)
		@tensor tmp[1 2; 4] := psi[i-1][1,2,3] * l[3,4] 
		psi[i-1] = tmp
	end
	_renormalize!(psi, psi[1], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

function _rightorth!(psi::FockMPS, alg::SVD, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	L = length(psi)
	maxerr = 0.
	for i in L:-1:2
		u, s, v, err = tsvd!(psi[i], (1,), (2, 3), trunc=trunc)
		psi[i] = v
		nr = _renormalize!(psi, s, normalize)
		rerror = sqrt(err * err / (nr * nr + err * err))
		(verbosity > 1) && println("SVD truncerror at bond $(i): ", rerror)
		u2 = u * Diagonal(s)
		@tensor tmp[-1 -2; -3] := psi[i-1][-1, -2, 1] * u2[1, -3]
		psi[i-1] = tmp
		psi.s[i] = s
		maxerr = max(maxerr, err)
	end
	(verbosity > 0) && println("Max SVD truncerror in rightorth: ", maxerr)
	_renormalize!(psi, psi[1], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

function canonicalize!(psi::FockMPS; alg::Orthogonalize = Orthogonalize(trunc=DefaultTruncation, normalize=false))
	alg.normalize && @warn "canonicalize with renormalization not recommanded for FockMPS"
	L = length(psi)
	_leftorth!(psi, QR(), NoTruncation(), alg.normalize, alg.verbosity)
	_rightorth!(psi, alg.orth, alg.trunc, alg.normalize, alg.verbosity)
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


function _renormalize!(psi, r, normalize::Bool)
	nr = norm(r)
	if nr != zero(nr)
		if !normalize
			_rescaling!(psi, nr)
		end
		lmul!(1/nr, r)  
  	end
  	return nr
end

function _renormalize_coeff!(psi, normalize::Bool)
	normalize && setscaling!(psi, 1)
end