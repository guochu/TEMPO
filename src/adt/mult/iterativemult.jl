# DMRG1 / DMRGAlgorithm 定义于 src/algorithms.jl。
#
# ALS 引擎（HadamardCache、iterative_compute!、leftsweep!/rightsweep!）来自
# FiniteMPSAlgorithms：ADT×ADT 的乘法物理指标共享，对应其 hadamard（逐点乘积）
# 问题，环境传递 ⟨z| x ⊙ y ⟩ 与旧版 updatemultleft/right 完全一致。
#
# TEMPO 特有行为保留在本 wrapper 中：
# * 初始猜测 `alg.initguess`（`:svd` 流式 SVD / `:pre` / `:rand`）；
# * 收敛判据采用 FiniteMPSAlgorithms 的 `iterative_compute!`
#   （相邻 sweep 末位损失的相对变化）；
# * finalize：一次 QR leftsweep + 带截断的末次 right sweep（旧版
#   `rightsweep_final!`，同时把归一化键谱写入 `z.s`）。

function iterativemult(x::ADT, y::ADT, alg::DMRG1)
    if alg.initguess == :svd
        z = _svd_guess(x, y, alg.trunc.D)
    elseif alg.initguess == :rand
        z = randomadt(promote_type(scalartype(x), scalartype(y)), phydims(x), D=alg.trunc.D)
        canonicalize!(z, alg=Orthogonalize(normalize=true))
    elseif alg.initguess == :pre
        z = changebond!(copy(x), alg.trunc.D)
        setscaling!(z, 1)
    else
        error("unsupported initguess $(alg.initguess)")
    end
    fmaalg = _fmadmrg1(alg)
    vz = z.parent
    cache = HadamardCache(x.parent, y.parent, vz)
    iterative_compute!(cache, fmaalg)
    _finalize!(cache, fmaalg, alg.trunc)
    setscaling!(z, scaling(x) * scaling(y))
    _rescaling!(z)
    return z
end

# finalize sweep：与旧版 rightsweep_final! 一致——QR 左扫重建左环境后，
# 从右到左以 `trunc` 做截断 SVD，并把归一化的键谱写入 Schmidt 值。
function _finalize!(m::HadamardCache, alg::FiniteMPSAlgorithms.DMRG1, trunc::TruncationScheme)
    L = length(m.bra)
    leftsweep!(m, alg)
    for site in L:-1:2
        mpsj = _reduce_hadamard_site(m.ketx[site], m.kety[site], m.hstorage[site], m.hstorage[site+1])
        u, s, v = tsvd!(mpsj, (1,), (2, 3), trunc=trunc)
        m.bra[site] = permute(v, (1, 2), (3,))
        if site == 2
            m.bra[1] = _contract_last(m.bra[1], u * Diagonal(s))
        end
        m.bra.s[site] = normalize!(s)
        m.hstorage[site] = _updateright(m.hstorage[site+1], m.bra[site], m.ketx[site], m.kety[site])
    end
    return m
end


# provide the initial guess（流式 SVD：逐 site 张量积经 truncdim(D) 截断）
_svd_guess(x::ADT, y::ADT, D::Int) = _svd_guess!(copy(x), y, D)
function _svd_guess!(x::ADT, y::ADT, D::Int)
    (length(x) == length(y)) || throw(DimensionMismatch())
     T = promote_type(scalartype(x), scalartype(y))
    left = ones(T, 1, 1, 1)
    tmp5 = n_fuse(_mult_site_n(x[1], y[1]), 3)
    @tensor tmp4[1,4;5,6] := left[1,2,3] * tmp5[2,3,4,5,6]
    trunc = truncdim(D)
    for i in 1:length(x)-1
        u, s, v = tsvd!(tmp4, (1,2), (3,4), trunc=trunc)
        x[i] = u
        _renormalize!(x, s, false)
        s2 = Matrix(Diagonal(s))
        @tensor r[1,3,4] := s2[1,2] * v[2,3,4]
        @tensor tmp1[1,5,4;2] := r[1,2,3] * y[i+1][3,4,5]
        @tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * x[i+1][4,5,6]
        tmp4 = n_fuse(tmp2, 2)

    end
    @tensor tmp[1,2;5] := tmp4[1,2,3,4] * conj(left[5,3,4])
    x[end] = tmp
    _rightorth!(x, SVD(), trunc, false, 0)
    setscaling!(x, 1)
    return x
end
