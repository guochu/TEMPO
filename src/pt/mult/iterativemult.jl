# DMRG1 / DMRGAlgorithm 定义于 src/algorithms.jl。
#
# ALS 引擎（MultCache、iterative_compute!、leftsweep!/rightsweep!）来自
# FiniteMPSAlgorithms：ProcessTensor×ProcessTensor 对应其 MPO×MPO 的 mult
# 问题，三链环境 ⟨ompo| mpo |impo⟩ 与旧版 updateleft/updateright 完全一致。
#
# TEMPO 特有行为保留在本 wrapper 中：
# * 初始猜测 `alg.initguess`（`:svd` 经 FMA 的 `svdguess_mult` 流式 SVD /
#   `:pre` / `:rand`）；
# * 收敛判据采用 FiniteMPSAlgorithms 的 `iterative_compute!`；
# * finalize：QR 左扫 + 末次 right sweep（共享的 `_finalize!`，见
#   fmabackend.jl；与 ADT 侧不同，PT 侧以 `NoTruncation()` 调用——保持旧行
#   为——只重新规范化并把归一化键谱写入 `z.s`）。

function iterativemult(x::ProcessTensor, y::ProcessTensor, alg::DMRG1)
    if alg.initguess == :svd
        z = ProcessTensor(svdguess_mult(x.parent, y.parent, alg.trunc.D))
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
    _finalize!(cache, fmaalg, NoTruncation())
    setscaling!(z, scaling(x) * scaling(y))
    _rescaling!(z)
    return z
end
