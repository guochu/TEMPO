

struct MPOMPOIterativeMultCache{_MPO, _IMPO, _OMPO, _H}
    ompo::_OMPO
	mpo::_MPO
	impo::_IMPO
	hstorage::_H
end


function iterativemult(x::ProcessTensor, y::ProcessTensor, alg::DMRGMultAlgorithm)
    # T = promote_type(eltype(mpo), eltype(mps))
    # mpsout = randommpo(T, ophysical_dimensions(mpo), iphysical_dimensions(mps), D=alg.D)
    # rightorth!(mpsout, alg=Orthogonalize(normalize=true))
    if alg.initguess == :svd
        z = _svd_guess(x, y, alg.D)
    elseif alg.initguess == :rand
        z = randompt(promote_type(scalartype(x), scalartype(y)), phydims(x), D=alg.D)
    elseif alg.initguess == :pre
        z = increase_bond!(copy(x), alg.D)
    else
        error("unsupported initguess $(alg.initguess)")
    end
    cache = mult_cache(z, x, y)
    deltas = compute!(cache, alg)
    z = cache.ompo
    setscaling!(z, scaling(x) * scaling(y))
    _rescaling!(z)
    return z
end

function mult_cache(z::ProcessTensor, x::ProcessTensor, y::ProcessTensor)
    @assert length(z) == length(x) == length(y)
    for i in 1:length(z)
        (phydim(z, i) == phydim(x, i) == phydim(y, i)) || throw(DimensionMismatch("phydim mismatch"))
    end
    # initialize Hstorage
    hstorage = init_hstorage_right(z, x, y)
    return MPOMPOIterativeMultCache(z, x, y, hstorage)
end

function finalize!(m::MPOMPOIterativeMultCache, alg::DMRGMultAlgorithm) end
function finalize!(m::MPOMPOIterativeMultCache, alg::DMRGMult1)
    leftsweep!(m, alg)
    rightsweep_final!(m, alg)
end

compute!(env::MPOMPOIterativeMultCache, alg::DMRGMultAlgorithm) = iterative_compute!(env, alg)

sweep!(m::MPOMPOIterativeMultCache, alg::DMRGMultAlgorithm) = vcat(leftsweep!(m, alg), rightsweep!(m, alg))



function leftsweep!(m::MPOMPOIterativeMultCache, alg::DMRGMult1)
    mpoA = m.impo
    mpo = m.mpo
    mpoB = m.ompo
    Cstorage = m.hstorage
    L = length(mpo)
    kvals = Float64[]
    workspace = scalartype(mpoB)[]
    for site in 1:L-1
        (alg.verbosity > 3) && println("Sweeping from left to right at bond: $site")
        mpsj = reduceH_single_site(mpoA[site], mpo[site], Cstorage[site], Cstorage[site+1])
        push!(kvals, norm(mpsj))
        (alg.verbosity > 1) && println("residual is $(kvals[end])...")
		q, r = tqr!(mpsj, (1,2,4), (3,), workspace)
        mpoB[site] = permute(q, (1,2,4,3))
        Cstorage[site+1] = updateleft(Cstorage[site], mpoB[site], mpo[site], mpoA[site])
    end
    return kvals	
end

function rightsweep!(m::MPOMPOIterativeMultCache, alg::DMRGMult1)
    mpoA = m.impo
    mpo = m.mpo
    mpoB = m.ompo
    Cstorage = m.hstorage
    L = length(mpo)
    kvals = Float64[]
    r = zeros(scalartype(mpoB), 0, 0)
    workspace = scalartype(mpoB)[]
    for site in L:-1:2
        (alg.verbosity > 3) && println("Sweeping from right to left at bond: $site.")
        mpsj = reduceH_single_site(mpoA[site], mpo[site], Cstorage[site], Cstorage[site+1])
        push!(kvals, norm(mpsj))
        (alg.verbosity > 1) && println("residual is $(kvals[end])...")

        r, mpoB[site] = tlq!(mpsj, (1,), (2,3,4), workspace)

        Cstorage[site] = updateright(Cstorage[site+1], mpoB[site], mpo[site], mpoA[site])
    end
    # println("norm of r is $(norm(r))")
    mpoB[1] = @tensor tmp[1,2,5,4] := mpoB[1][1,2,3,4] * r[3,5]
    return kvals	
end

function rightsweep_final!(m::MPOMPOIterativeMultCache, alg::DMRGMult1)
    mpoA = m.impo
    mpo = m.mpo
    mpoB = m.ompo
    Cstorage = m.hstorage
    L = length(mpo)
    kvals = Float64[]
    r = zeros(scalartype(mpoB), 0, 0)
    workspace = scalartype(mpoB)[]
    for site in L:-1:2
        (alg.verbosity > 3) && println("Sweeping from right to left at bond: $site.")
        mpsj = reduceH_single_site(mpoA[site], mpo[site], Cstorage[site], Cstorage[site+1])
        push!(kvals, norm(mpsj))
        (alg.verbosity > 1) && println("residual is $(kvals[end])...")

        r, s, mpoB[site] = tsvd!(mpsj, (1,), (2,3,4), workspace)

        if site == 2
            r = r * Diagonal(s)
            mpoB[1] = @tensor tmp[1,2,5,4] := mpoB[1][1,2,3,4] * r[3,5]
        end
        mpoB.s[site] = normalize!(s)


        Cstorage[site] = updateright(Cstorage[site+1], mpoB[site], mpo[site], mpoA[site])
    end
    # println("norm of r is $(norm(r))")
    # mpoB[1] = @tensor tmp[1,2,5,4] := mpoB[1][1,2,3,4] * r[3,5]
    return kvals    
end

_svd_guess(x::ProcessTensor, y::ProcessTensor, D::Int) = _svd_guess!(copy(x), y, D)
function _svd_guess!(x::ProcessTensor, y::ProcessTensor, D::Int)
    (length(x) == length(y)) || throw(DimensionMismatch())
    T = promote_type(scalartype(x), scalartype(y))
    left = ones(T, 1, 1, 1)
    tmp5 = _mult_mpo_sitetensor(x[1], y[1])
    @tensor tmp4[1,4,5,6,7] := left[1,2,3] * tmp5[2,3,4,5,6,7]
    trunc = truncdim(D)
    for i in 1:length(x)-1
        u, s, v = tsvd!(tmp4, (1,2,5), (3,4), trunc=trunc)
        x[i] = permute(u, (1,2,4,3))
        _renormalize!(x, s, false)
        s2 = Matrix(Diagonal(s))
        @tensor r[1,3,4] := s2[1,2] * v[2,3,4]
        @tensor tmp1[1,6,5,4,2] := r[1,2,3] * y[i+1][3,4,5,6]
        @tensor tmp4[1,6,7,3,2] := tmp1[1,2,3,4,5] * x[i+1][5,6,7,4]
    end
    @tensor tmp[1,2,6,5] := tmp4[1,2,3,4,5] * left[6,3,4]
    x[end] = tmp
    _rightorth!(x, SVD(), trunc, false, 0)
    setscaling!(x, 1)
    return x
end

function reduceH_single_site(A::DenseMPOTensor, m::DenseMPOTensor, cleft::DenseMPSTensor, cright::DenseMPSTensor)
	@tensor tmp[1,7,9,6] := ((cleft[1,2,3] * A[3,4,5,6]) * m[2,7,8,4]) * cright[9,8,5]
    return tmp
end


function init_hstorage_right(B::ProcessTensor, mpo::ProcessTensor, A::ProcessTensor)
    @assert length(B) == length(mpo) == length(A)
    L = length(mpo)
    T = scalartype(B)
    hstorage = Vector{Array{T, 3}}(undef, L+1)
    hstorage[1] = ones(1,1,1)
    hstorage[L+1] = ones(1,1,1)
    for i in L:-1:2
        hstorage[i] = updateright(hstorage[i+1], B[i], mpo[i], A[i])
    end
    return hstorage
end


function updateleft(cleft::DenseMPSTensor, B::DenseMPOTensor, m::DenseMPOTensor, A::DenseMPOTensor)
    @tensor tmp[9,8,5] := ((cleft[1,2,3] * A[3,4,5,6]) * m[2,7,8,4]) * conj(B[1,7,9,6])
    return tmp
end

function updateright(cright::DenseMPSTensor, B::DenseMPOTensor, m::DenseMPOTensor, A::DenseMPOTensor)
    @tensor tmp[1,7,9] := ((conj(B[1,2,3,4]) * cright[3,5,6]) * m[7,2,5,8] ) * A[9,8,6,4]
    return tmp
end

# function _mult_site_n(xj::DenseMPOTensor, yj::DenseMPOTensor)
#     @tensor r[1,5,2,3,6,7] := xj[1,2,3,4] * yj[5,4,6,7]
#     return r
# end
