# DMRG1 / DMRGAlgorithm 定义于 src/algorithms.jl。
#
# ALS 引擎（HadamardCache、iterative_compute!、leftsweep!/rightsweep!）来自
# FiniteMPSAlgorithms：ADT×ADT 的乘法物理指标共享，对应其 hadamard（逐点乘积）
# 问题，环境传递 ⟨z| x ⊙ y ⟩ 与旧版 updatemultleft/right 完全一致。
#
# TEMPO 特有行为保留在本 wrapper 中：
# * 初始猜测 `alg.initguess`（`:svd` 经 FMA 的 `svdguess_hadamard` 流式 SVD /
#   `:pre` / `:rand`）；
# * 收敛判据采用 FiniteMPSAlgorithms 的 `iterative_compute!`；
# * finalize：对输出链做 FMA 的 `canonicalize!`（QR 左扫 + 以 `alg.trunc`
#   截断的 SVD 右扫，见 fmabackend.jl），键谱写入 `z.s`。

function iterativemult(x::ADT, y::ADT, alg::DMRG1)
    (length(x) == length(y)) || throw(DimensionMismatch())
    for i in 1:length(x)
        (phydim(x, i) == phydim(y, i)) || throw(DimensionMismatch("phydim mismatch"))
    end
    if alg.initguess == :svd
        z = ADT(svdguess_hadamard(x.parent, y.parent, alg.trunc.D))
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
    cache = HadamardCache(x.parent, y.parent, z.parent)
    iterative_compute!(cache, fmaalg)
    _finalize!(cache, fmaalg, alg.trunc)
    # `_finalize!`（canonicalize!）把链的范数因子折叠进 `scaling(z)`，与
    # svdmult 的 `setscaling!(x, scaling(x) * scaling(y))` 同构：输出的绝对
    # 幅值 = (scaling(z)·scaling(x)·scaling(y))^L × 单位规范链。
    setscaling!(z, scaling(z) * scaling(x) * scaling(y))
    _rescaling!(z)
    return z
end
