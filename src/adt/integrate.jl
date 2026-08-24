
"""
	integrate(x::ADT)

Fully contract the state represented by the MPS and return the resulting overall scalar value (equivalent to the total coefficient obtained by summing over all sites).

# Arguments
- `x::ADT`: MPS to contract

# Returns
Scalar: result of contracting the whole MPS.
"""
function integrate(x::ADT)
	L = length(x)
	sca = scaling(x)
	v = dropdims(sum(x[L], dims=2), dims=(2,3)) * sca
	for i in L-1:-1:1
		tmp = dropdims(sum(x[i], dims=2), dims=2) * sca
		v = tmp * v
	end
	return only(v)
end


"""
	integrate(x::ADT, y::ADT)

Compute the overlap scalar between two MPS of equal length.

# Arguments
- `x::ADT`: first MPS
- `y::ADT`: second MPS, must have the same length as `x`

# Returns
Scalar: the inner product obtained by contracting `x` with `y`.
"""
function integrate(x::ADT, y::ADT)
	(length(x) == length(y)) || throw(DimensionMismatch("adt size mismatch"))
	sca = scaling(x) * scaling(y)
	L = length(x)
	@tensor v[1,4] := sca * x[L][1,2,3] * y[L][4,2,3]
	for i in L-1:-1:1
		@tensor tmp[1,4] := sca * x[i][1,2,3] * y[i][4,2,5] * v[3,5]
		v = tmp
	end
	return tr(v)
end