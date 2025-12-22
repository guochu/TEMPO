function integrate(lat::ImagPTLattice, x::ProcessTensor)
	(length(lat) == length(x)) || throw(DimensionMismatch("lattice size mismatch with PT size"))
    L = length(x)
    sca = scaling(x)
    v = dropdims(x[L], dims=3) * sca
    for i in L-1:-1:1
        @tensor tmp[1,5,4] := sca * x[i][1,2,3,4] * v[3,5,2]
        v = tmp
    end
    v = dropdims(v, dims=1)
    return tr(v)
end

function integrate(lat::RealPTLattice, x::ProcessTensor)
	(length(lat) == length(x)) || throw(DimensionMismatch("lattice size mismatch with PT size"))
    L = lat.N
    sca2 = scaling(x)^2
    pos1, pos2 = index(lat, L, branch=:+), index(lat, L, branch=:-)
    @assert pos1+1 == pos2 == 2
    v = dropdims(x[pos1], dims=1)
    @tensor tmp[3,5,4] := v[1,2,3] * x[pos2][2,1,4,5] * sca2
    for i in L-1:-1:1
    	pos1, pos2 = index(lat, i, branch=:+), index(lat, i, branch=:-)
    	@assert pos1+1 == pos2
    	@tensor v[5,7,6] := tmp[1,2,3] * x[pos1][3,1,4,5] * x[pos2][4,2,6,7] * sca2
    	tmp = v
    end
    return tr(dropdims(tmp, dims=3))
end

function integrate(lat::MixedPTLattice, x::ProcessTensor)
	(length(lat) == length(x)) || throw(DimensionMismatch("lattice size mismatch with PT size"))
	sca = scaling(x)
	sca2 = sca^2
	pos = index(lat, lat.Nτ, branch=:τ)
	@assert pos == 1
	v = dropdims(x[pos], dims=1) * sca
	for i in 2:lat.Nτ
		pos = index(lat, i, branch=:τ)
		@tensor tmp[1,4,5] := sca * v[1,2,3] * x[i][2,3,4,5]
		v = tmp
	end
	pos1, pos2 = index(lat, 1, branch=:-), index(lat, 1, branch=:+)
	@assert pos1 + 1 == pos2
	@tensor v2[4,6,7] := v[1,2,3] * x[pos1][2,4,5,3] * x[pos2][5,6,7,1] * sca2
	for i in 2:lat.Nt
		pos1, pos2 = index(lat, i, branch=:-), index(lat, i, branch=:+)
		@tensor tmp[4,6,7] := v2[1,2,3] * x[pos1][3,4,5,1] * x[pos2][5,6,7,2] * sca2
		v2 = tmp
	end
	return tr(dropdims(v2, dims=3))
end


# function ADT(x::ProcessTensor; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
#     L = length(x)
#     (L > 1) || throw(DimensionMismatch("ProcessTensor size should be larger than 1"))
#     T = scalartype(x)
#     data = Vector{Array{T, 3}}(undef, L+1)
#     d = size(x[1], 2)
#     c = copy_tensor3(T, d)
#     workspace = T[]
#     q, r = tqr!(x[1], (1,4), (2, 3), workspace)
#     data[1] = q
#     for i in 2:L
#         @tensor tmp2[1,7,4,5] := r[1,2,3] * x[i][3,4,5,6] * c[7,2,6]
#         q, r = tqr!(tmp2, (1,2), (3,4), workspace)
#         data[i] = q
#     end
#     data[L+1] = r

#     y = ADT(data)
#     _rightorth!(y, SVD(), trunc, false, verbosity)
#     setscaling!(y, scaling(x)^(L/L+1) * scaling(y))
#     return y
# end


# function copy_tensor3(::Type{T}, d::Int) where {T<:Number}
#     m = zeros(T, d, d, d)
#     for i in 1:d
#         m[i,i,i] = 1
#     end
#     return m
# end