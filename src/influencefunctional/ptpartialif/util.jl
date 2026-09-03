"""
	partialif_naive(lattice::AbstractPTLattice, rowind::ContourIndex, corr::AbstractCorrelationFunction, hyb::NonAdditiveHyb; trunc::TruncationScheme=DefaultITruncation)

Naive partial influence functional construction on PT lattices: apply the exp(ηᵢⱼ op⊗op) gates successively for each column index (`_get_contour_op` handles the transpose of op according to the contour branch), orthogonalizing/truncating along the way.

# Arguments
- `lattice::AbstractPTLattice`: PT contour lattice.
- `rowind::ContourIndex`: contour index (row) associated with the partial IF.
- `corr::AbstractCorrelationFunction`: bath correlation function.
- `hyb::NonAdditiveHyb`: non-additive (symmetric) system-bath coupling.
- `trunc::TruncationScheme=DefaultITruncation`: truncation scheme used during orthogonalization.

# Returns
The partial IF, represented as a `ProcessTensor`.
"""
function partialif_naive(lattice::AbstractPTLattice, rowind::ContourIndex, corr::AbstractCorrelationFunction, hyb::NonAdditiveHyb; 
							trunc::TruncationScheme=DefaultITruncation)
	z = hyb.op
	(lattice.d == size(z, 1) == size(z, 2)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))
	# d = lattice.d
	# z2 = z * z
	# zz = kron(z, z)

	b1 = branch(rowind)
	i = rowind.j
	# pos1 = index(lattice, i, branch=b1)
		
	T = promote_type(scalartype(lattice), scalartype(hyb), scalartype(corr))
	tmp = vacuumstate(T, lattice)
	orth = Orthogonalize(SVD(), trunc)
	for b2 in branches(lattice)
		k2 = (b2 == :τ) ? lattice.Nτ : lattice.Nt
		for j in 1:k2
			coef = index(corr, i, j, b1=b1, b2=b2)
			ind2 = ContourIndex(j, b2) 
			# pos2 = lattice[ind2]
			# if pos1 == pos2
			# 	m = exp(coef .* z2)
			# 	# t = ContourOperator(ind1, m)
			# 	t = FockTermS(pos1, m)
			# else
			# 	m = exp(coef .* zz)
			# 	t = FockTermS((pos1, pos2), reshape(m, (d,d,d,d)))
			# end
			t = _get_contour_op(lattice, rowind, ind2, z, coef)
			apply!(t, tmp)
			canonicalize!(tmp, alg=orth)
		end
	end
	return tmp
end

function _get_contour_op(lattice, ind1::ContourIndex, ind2::ContourIndex, z::AbstractMatrix, coef)
	d = lattice.d
	z1 = z
	z2 = z
	if branch(ind1) == :-
		z1 = transpose(z)
	end
	if branch(ind2) == :-
		z2 = transpose(z)
	end
	pos1, pos2 = lattice[ind1], lattice[ind2]
	if pos1 == pos2
		m = exp(coef .* z1 * z2) 
		t = FockTermS(pos1, m)
	else
		zz = kron(z2, z1)
		m = exp(coef .* zz)
		t = FockTermS((pos1, pos2), reshape(m, (d,d,d,d)))
	end
	return t
end

"""
	partialif(lattice::AbstractPTLattice, rowind::ContourIndex, corr::AbstractCorrelationFunction, hyb::NonAdditiveHyb)

Construct the partial influence functional (partial IF) associated with the contour index `rowind` as a `ProcessTensor` with exact bond dimension `d` (the physical dimension), generalizing the algorithm of Strathearn et al. (2018) to non-additive couplings.

For a Hermitian coupling operator z = V·Diagonal(λ)·V' (zᵀ = conj(V)·Diagonal(λ)·conj(V)'), the contour gates are exp(ηᵢⱼ·op(b₁)@row ⊗ op(b₂)@col) with op(:+) = op(:τ) = z and op(:-) = zᵀ. All row-site gates carry the same operator op(b₁), so the row operators collapse onto its spectral projectors, and the partial IF is

	P = Σₐ e^{ηᵢᵢλₐ²}P₍b₁₎⁽ᵃ⁾@row · ∏ⱼ exp(ηᵢⱼλₐ·op(b₂))@colⱼ,

i.e. an MPO whose bond index is the eigenindex `a`.

# Arguments
- `lattice::AbstractPTLattice`: PT contour lattice.
- `rowind::ContourIndex`: contour index (row) associated with the partial IF.
- `corr::AbstractCorrelationFunction`: bath correlation function.
- `hyb::NonAdditiveHyb`: non-additive Hermitian system-bath coupling.

# Returns
The partial IF, represented as a `ProcessTensor`.
"""
function partialif(lattice::AbstractPTLattice, rowind::ContourIndex, corr::AbstractCorrelationFunction, hyb::NonAdditiveHyb)
	z = hyb.op
	(lattice.d == size(z, 1) == size(z, 2)) || throw(DimensionMismatch("lattice.d mismatch with hyb.d"))

	b1 = branch(rowind)
	i = rowind.j
	pos1 = index(lattice, i, branch=b1)

	T = promote_type(scalartype(lattice), scalartype(hyb), scalartype(corr))
	colsites = Tuple{Int, Symbol}[]
	coefs = T[]
	for b2 in branches(lattice)
		k2 = (b2 == :τ) ? lattice.Nτ : lattice.Nt
		for j in 1:k2
			pos2 = index(lattice, j, branch=b2)
			coef = index(corr, i, j, b1=b1, b2=b2)
			push!(colsites, (pos2, b2))
			push!(coefs, coef)
		end
	end
	return partialif_densempo(length(lattice), lattice.d, pos1, colsites, coefs, z, b1)
end

# generalization of the StrathearnLovett2018 algorithm to matrix (non-additive) couplings
function partialif_densempo(L::Int, d::Int, row::Int, colsites::Vector{Tuple{Int, Symbol}},
							coefs::Vector{<:Number}, op::AbstractMatrix{<:Number}, rowbranch::Symbol)
	@assert length(colsites) == length(coefs) == L > 1
	p = sortperm(first.(colsites))
	colsites = colsites[p]
	coefs = coefs[p]
	T = promote_type(scalartype(op), scalartype(coefs))

	# op is Hermitian by construction of NonAdditiveHyb: op = V·Λ·V' and opᵀ = conj(V)·Λ·conj(V)'
	F = (eltype(op) <: Real) ? eigen(Symmetric(op)) : eigen(Hermitian(op))
	V = F.vectors
	Vc = conj(V)
	λ = Vector{T}(F.values)

	# all row-site gates carry the same operator op(rowbranch), so the row operators
	# collapse onto its spectral projectors and the bond index is the eigenindex
	Wr = (rowbranch == :-) ? Vc : V
	rowcontent(a) = exp(coefs[findfirst(x->x[1] == row, colsites)] * λ[a]^2) .*
					(Wr[:, a:a] * adjoint(Wr)[a:a, :])
	# column site j carries op(b₂): exp(ηᵢⱼλₐ·op(b₂))
	function colcontent(s, a)
		Wc = (colsites[s][2] == :-) ? Vc : V
		return Wc * Diagonal(exp.(coefs[s] * λ[a] .* λ)) * adjoint(Wc)
	end
	content(s, a) = (s == findfirst(x->x[1] == row, colsites)) ? rowcontent(a) : colcontent(s, a)

	mpsdata = Vector{Array{T, 4}}(undef, L)
	for s in 1:L
		if s == 1
			tmp = zeros(T, 1, d, d, d)
			for a in 1:d
				tmp[1, :, a, :] .= content(s, a)
			end
		elseif s == L
			tmp = zeros(T, d, d, 1, d)
			for a in 1:d
				tmp[a, :, 1, :] .= content(s, a)
			end
		else
			tmp = zeros(T, d, d, d, d)
			for a in 1:d
				tmp[a, :, a, :] .= content(s, a)
			end
		end
		mpsdata[s] = tmp
	end
	return ProcessTensor(mpsdata)
end


