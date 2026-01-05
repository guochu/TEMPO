# function integrate(lat::ImagPTLattice, x::ProcessTensor)
# 	(length(lat) == length(x)) || throw(DimensionMismatch("lattice size mismatch with PT size"))
#     L = length(x)
#     sca = scaling(x)
#     v = dropdims(x[L], dims=3) * sca
#     for i in L-1:-1:1
#         @tensor tmp[1,5,4] := sca * x[i][1,2,3,4] * v[3,5,2]
#         v = tmp
#     end
#     v = dropdims(v, dims=1)
#     return tr(v)
# end

function integrate(lat::ImagPTLattice, x::ProcessTensor)
	(length(lat) == length(x)) || throw(DimensionMismatch("lattice size mismatch with PT size"))
    L = length(x)
    sca = scaling(x)
    v = dropdims(x[L], dims=3) * sca
    for i in L-1:-1:1
        @tensor tmp[1,2,5] := sca * x[i][1,2,3,4] * v[3,4,5]
        v = tmp
    end
    v = dropdims(v, dims=1)
    return tr(v)
end

# zipup algorithm
function integrate(lat::ImagPTLattice, x::ProcessTensor, y::ProcessTensor)
	(length(lat) == length(x) == length(y)) || throw(DimensionMismatch("lattice size mismatch with PT size"))
	L = length(x)
	sca = scaling(x) * scaling(y)
	v1 = dropdims(x[L], dims=3)
	v2 = dropdims(y[L], dims=3)
	@tensor v[1,4,2,5] := v1[1,2,3] * v2[4,3,5] * sca
	for i in L-1:-1:1
		@tensor tmp[7,1,8,6] := sca * y[i][1,2,3,4] * v[5,3,4,6] * x[i][7,8,5,2]
		v = tmp
	end
	v = dropdims(v, dims=(1,2))
	return tr(v)
end

# function integrate(lat::RealPTLattice, x::ProcessTensor)
# 	(length(lat) == length(x)) || throw(DimensionMismatch("lattice size mismatch with PT size"))
#     L = lat.N
#     sca2 = scaling(x)^2
#     pos1, pos2 = index(lat, L, branch=:+), index(lat, L, branch=:-)
#     @assert pos1+1 == pos2 == 2
#     v = dropdims(x[pos1], dims=1)
#     @tensor tmp[3,5,4] := v[1,2,3] * x[pos2][2,1,4,5] * sca2
#     for i in L-1:-1:1
#     	pos1, pos2 = index(lat, i, branch=:+), index(lat, i, branch=:-)
#     	@assert pos1+1 == pos2
#     	@tensor v[5,7,6] := tmp[1,2,3] * x[pos1][3,1,4,5] * x[pos2][4,2,6,7] * sca2
#     	tmp = v
#     end
#     return tr(dropdims(tmp, dims=3))
# end
integrate(lat::RealPTLattice, x::ProcessTensor) = tr(rdm(lat, x))
integrate(lat::RealPTLattice, x::ProcessTensor, y::ProcessTensor) = tr(rdm(lat, x, y))


"""
	rdm(lat::RealPTLattice, x::ProcessTensor)

The final output quantum state
"""
function rdm(lat::RealPTLattice, x::ProcessTensor)
	(length(lat) == length(x)) || throw(DimensionMismatch("lattice size mismatch with PT size"))
    L = lat.N
    sca2 = scaling(x)^2
    pos1, pos2 = index(lat, 1, branch=:+), index(lat, 1, branch=:-)
    @assert pos1+1 == pos2 == length(lat)
    v = dropdims(x[pos2], dims=3)
    @tensor tmp[1,2,5] := x[pos1][1,2,3,4] * v[3,5,4] * sca2
    for i in 2:L
    	pos1, pos2 = index(lat, i, branch=:+), index(lat, i, branch=:-)
    	@assert pos1+1 == pos2
    	@tensor v[1,2,5] := sca2 * x[pos1][1,2,3,4] * x[pos2][3,5,6,7] * tmp[6,4,7] 
    	tmp = v
    end
    return dropdims(tmp, dims=1)
end

function rdm(lat::RealPTLattice, x::ProcessTensor, y::ProcessTensor)
	(length(lat) == length(x) == length(y)) || throw(DimensionMismatch("lattice size mismatch with PT size"))
    L = lat.N
    sca2 = (scaling(x) * scaling(y))^2
    pos1, pos2 = index(lat, 1, branch=:+), index(lat, 1, branch=:-)
    @assert pos1+1 == pos2 == length(lat)
    v1 = dropdims(x[pos2], dims=3)
    v2 = dropdims(y[pos2], dims=3)
    @tensor v[1,4,2,5] := v1[1,2,3] * v2[4,3,5] 
    @tensor tmp[7,1,8,6] := y[pos1][1,2,3,4] * v[5,3,6,4] * x[pos1][7,8,5,2] * sca2
    for i in 2:L
    	pos1, pos2 = index(lat, i, branch=:+), index(lat, i, branch=:-)
    	@assert pos1+1 == pos2
    	@tensor v[7,1,8,6] := y[pos2][1,2,3,4] * tmp[5,3,6,4] * x[pos2][7,8,5,2]
    	@tensor tmp[7,1,8,6] := y[pos1][1,2,3,4] * v[5,3,6,4] * x[pos1][7,8,5,2] * sca2
    end
    return dropdims(tmp, dims=(1,2))
end

function quantummap(lat::RealPTLattice, x::ProcessTensor)
	(length(lat) == length(x)) || throw(DimensionMismatch("lattice size mismatch with PT size"))
    L = lat.N
    sca2 = scaling(x)^2
    pos1, pos2 = index(lat, 1, branch=:+), index(lat, 1, branch=:-)
    @assert pos1+1 == pos2 == length(lat)
    v = dropdims(x[pos2], dims=3)
    @tensor tmp[1,2,5,4,6] := x[pos1][1,2,3,4] * v[3,5,6] * sca2
    for i in 2:L
    	pos1, pos2 = index(lat, i, branch=:+), index(lat, i, branch=:-)
    	@assert pos1+1 == pos2
    	@tensor v[1,2,5,8,9] := sca2 * x[pos1][1,2,3,4] * x[pos2][3,5,6,7] * tmp[6,4,7,8,9] 
    	tmp = v
    end
    return dropdims(tmp, dims=1)
	
end

function integrate(lat::MixedPTLattice, x::ProcessTensor)
	(length(lat) == length(x)) || throw(DimensionMismatch("lattice size mismatch with PT size"))
	sca = scaling(x)
	sca2 = sca^2
	pos = index(lat, lat.Nτ, branch=:τ)
	@assert pos == 1
	v = dropdims(x[pos], dims=1) * sca
	for i in lat.Nτ-1:-1:1
		pos = index(lat, i, branch=:τ)
		@tensor tmp[1,4,5] := sca * v[1,2,3] * x[pos][2,3,4,5]
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

function integrate(lat::MixedPTLattice, x::ProcessTensor, y::ProcessTensor)
	(length(lat) == length(x) == length(y)) || throw(DimensionMismatch("lattice size mismatch with PT size"))
	sca = scaling(x) * scaling(y)
	sca2 = sca^2
	pos = index(lat, lat.Nτ, branch=:τ)
	@assert pos == 1
	v1 = dropdims(x[pos], dims=1) 
	v2 = dropdims(y[pos], dims=1)
	@tensor v[1,2,4,5] := v1[1,2,3] * v2[3,4,5] * sca
	for i in lat.Nτ-1:-1:1
		pos = index(lat, i, branch=:τ)
		@tensor tmp[1,5,7,8] := sca * v[1,2,3,4] * x[pos][2,4,5,6] * y[pos][3,6,7,8]
		v = tmp
	end
	pos1, pos2 = index(lat, 1, branch=:-), index(lat, 1, branch=:+)
	@assert pos1 + 1 == pos2
	@tensor tmp[1,7,8,6] := v[1,2,3,4] * y[pos1][3,5,6,4] * x[pos1][2,7,8,5]
	@tensor v2[2,7,8,6] := tmp[1,2,3,4] * y[pos2][4,5,6,1] * x[pos2][3,7,8,5] * sca2
	for i in 2:lat.Nt
		pos1, pos2 = index(lat, i, branch=:-), index(lat, i, branch=:+)
		@tensor tmp[2,7,8,6] := v2[1,2,3,4] * y[pos1][4,5,6,1] * x[pos1][3,7,8,5]
		@tensor v2[2,7,8,6] := tmp[1,2,3,4] * y[pos2][4,5,6,1] * x[pos2][3,7,8,5] * sca2
	end
	return tr(dropdims(v2, dims=(3,4)))
end
