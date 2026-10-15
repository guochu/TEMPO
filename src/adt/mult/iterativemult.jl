# DMRG1 / DMRG2 / DMRGAlgorithm 定义于 src/algorithms.jl。
#
# ALS 引擎（HadamardCache、iterative_compute!、leftsweep!/rightsweep!）来自
# FiniteMPSAlgorithms：ADT×ADT 的乘法物理指标共享，对应其 hadamard（逐点乘积）
# 问题，环境传递 ⟨z| x ⊙ y ⟩ 与旧版 updatemultleft/right 完全一致。
#
# TEMPO 特有行为保留在本 wrapper 中：
# * 初始猜测 `alg.initguess`（`:svd` 经 FMA 的 `svdguess_hadamard` 流式 SVD /
#   `:pre` / `:rand`；三者均输出右正则链——FMA 的 cache 构造要求右正则 bra
#   且不再代为规范）；DMRG2 的 `:svd` 同上，`:pre` 只复制键维较大的输入
#   （两站点扫描自行适配键维，不调 `changebond!` 升键维；初值须已右正则）；
# * 收敛判据采用 FiniteMPSAlgorithms 的 `iterative_compute!`（返回
#   `ALSConvergenceInfo`；未收敛时 FMA 已按 verbosity 自行告警）；
# * finalize：对输出链做 FMA 的 `truncate!`（QR 左扫 + 以 `alg.trunc`
#   截断的 SVD 右扫，见 fmabackend.jl），键谱写入 `z.s`。

function iterativemult(x::ADT, y::ADT, alg::DMRG1)
    (length(x) == length(y)) || throw(DimensionMismatch())
    for i in 1:length(x)
        (phydim(x, i) == phydim(y, i)) || throw(DimensionMismatch("phydim mismatch"))
    end
    if alg.initguess == :svd
        z = ADT(svdguess_hadamard(x.parent, y.parent, alg.trunc.D))
    elseif alg.initguess == :rand
        # randomadt is already right-canonical (FMA's randommps); no re-gauging needed
        z = randomadt(promote_type(scalartype(x), scalartype(y)), phydims(x), D=alg.trunc.D)
    elseif alg.initguess == :pre
        # 键维更大的输入包含的乘积结构信息更多，以其为初值收敛更快；键维相同时用 x
        z = changebond!(copy(bonddim(x) ≥ bonddim(y) ? x : y), alg.trunc.D)
        setscaling!(z, 1)
    else
        error("unsupported initguess $(alg.initguess)")
    end
    fmaalg = _fmadmrg1(alg)
    # HadamardCache(bra, ketx, kety)：bra = z，ketx = x，kety = y
    cache = HadamardCache(z.parent, x.parent, y.parent)
    info = iterative_compute!(cache, fmaalg)
    _finalize!(cache, alg)
    # `_finalize!`（truncate!）把链的范数因子折叠进 `scaling(z)`，与
    # svdmult 的 `setscaling!(x, scaling(x) * scaling(y))` 同构：输出的绝对
    # 幅值 = (scaling(z)·scaling(x)·scaling(y))^L × 单位规范链。
    setscaling!(z, scaling(z) * scaling(x) * scaling(y))
    # 首位张量范数折入 scaling（数值区间管理；与 FiniteMPSAlgorithms 的
    # `_renormalize!` 一致）
    _renormalize!(z, z[1], false)
    return z
end

# 两站点路线：`:svd` 初值取精确乘积的流式 SVD（键维上限 trunc.D，右正则链）；
# `:pre` 只复制键维较大的输入（键维由两站点扫描自适应，无需 changebond! 升维；
# 初值须已是右正则链，copy 保持其规范）。其余流程（scaling 记账、finalize）与
# DMRG1 路线一致。
function iterativemult(x::ADT, y::ADT, alg::DMRG2)
    (length(x) == length(y)) || throw(DimensionMismatch())
    for i in 1:length(x)
        (phydim(x, i) == phydim(y, i)) || throw(DimensionMismatch("phydim mismatch"))
    end
    if alg.initguess == :svd
        z = ADT(svdguess_hadamard(x.parent, y.parent, alg.trunc.D))
    elseif alg.initguess == :pre
        z = copy(bonddim(x) ≥ bonddim(y) ? x : y)
        setscaling!(z, 1)
    else
        error("unsupported initguess $(alg.initguess)")
    end
    fmaalg = _fmadmrg2(alg)
    cache = HadamardCache(z.parent, x.parent, y.parent)
    info = iterative_compute!(cache, fmaalg)
    _finalize!(cache, alg)
    setscaling!(z, scaling(z) * scaling(x) * scaling(y))
    _renormalize!(z, z[1], false)
    return z
end
