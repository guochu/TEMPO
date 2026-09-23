
"""
	integrate(x::ADT)

Fully contract the state represented by the MPS and return the resulting overall scalar value (equivalent to the total coefficient obtained by summing over all sites).

Delegates to FiniteMPSAlgorithms' `sum(::CanonicalMPS)`（所有物理指标与全一向量
收缩，含 `scaling^L` 因子）。

# Arguments
- `x::ADT`: MPS to contract

# Returns
Scalar: result of contracting the whole MPS.
"""
integrate(x::ADT) = sum(x.parent)


"""
	integrate(x::ADT, y::ADT)

Sum of the amplitudes of the pointwise (Hadamard) product `x ⊙ y` — equivalent to
`sum(x ⊙ y)`, realized lazily: the physical indices of `x` and `y` are shared (summed
over, no conjugation) while the bonds are contracted site by site.

# Arguments
- `x::ADT`: first MPS
- `y::ADT`: second MPS, must have the same length as `x`

# Returns
Scalar: the summed amplitude of the pointwise product.
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
