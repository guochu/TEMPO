# DMRG1 / DMRGAlgorithm 定义于 src/algorithms.jl。
#
# ALS 引擎（MultCache、iterative_compute!、leftsweep!/rightsweep!）来自
# FiniteMPSAlgorithms：ProcessTensor×ProcessTensor 对应其 MPO×MPO 的 mult
# 问题，三链环境 ⟨ompo| mpo |impo⟩ 与旧版 updateleft/updateright 完全一致。
#
# TEMPO 特有行为保留在本 wrapper 中：
# * 初始猜测 `alg.initguess`（`:svd` 流式 SVD / `:pre` / `:rand`）；
# * 收敛判据采用 FiniteMPSAlgorithms 的 `iterative_compute!`；
# * finalize：一次 QR leftsweep + 末次 right sweep（旧版 `rightsweep_final!`；
#   与 ADT 侧不同，PT 侧的末次 sweep 不带截断——保持旧行为——只重新规范
#   化并把归一化键谱写入 `z.s`）。

function iterativemult(x::ProcessTensor, y::ProcessTensor, alg::DMRG1)
    if alg.initguess == :svd
        z = _svd_guess(x, y, alg.trunc.D)
    elseif alg.initguess == :rand
        z = randompt(promote_type(scalartype(x), scalartype(y)), phydims(x), D=alg.trunc.D)
        canonicalize!(z, alg=Orthogonalize(normalize=true))
    elseif alg.initguess == :pre
        z = changebond!(copy(x), alg.trunc.D)
        setscaling!(z, 1)
    else
        error("unsupported initguess $(alg.initguess)")
    end
    fmaalg = _fmadmrg1(alg)
    vz = z.parent
    # MultCache(H, ket, bra)：H = mpo（x），ket = impo（y），bra = ompo（z）
    cache = MultCache(x.parent, y.parent, vz)
    iterative_compute!(cache, fmaalg)
    _finalize!(cache, fmaalg)
    setscaling!(z, scaling(x) * scaling(y))
    _rescaling!(z)
    return z
end

# finalize sweep：与旧版 rightsweep_final! 一致——QR 左扫重建左环境后，
# 从右到左做一次无截断的 SVD 重新规范化，并把归一化键谱写入 Schmidt 值。
function _finalize!(m::MultCache, alg::FiniteMPSAlgorithms.DMRG1)
    L = length(m.bra)
    leftsweep!(m, alg)
    for site in L:-1:2
        mpsj = _reduce_site(m.ket[site], m.H[site], m.hstorage[site], m.hstorage[site+1])
        u, s, v = tsvd!(mpsj, (1,), (2, 3, 4))
        m.bra[site] = v
        if site == 2
            m.bra[1] = _contract_last(m.bra[1], u * Diagonal(s))
        end
        m.bra.s[site] = normalize!(s)
        m.hstorage[site] = _env_updateright(m.hstorage[site+1], m.bra[site], m.H[site], m.ket[site])
    end
    return m
end


# provide the initial guess（流式 SVD：逐 site 张量积经 truncdim(D) 截断）
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
