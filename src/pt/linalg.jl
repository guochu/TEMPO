


# 精确（未压缩）算符链乘积：委托给 FiniteMPSAlgorithms 的 scaling-safe 乘法
# （CanonicalMPO × CanonicalMPO 的结果自带 scaling(x)·scaling(y)）
function Base.:*(x::ProcessTensor, y::ProcessTensor)
    @assert !isempty(x)
    (length(x) == length(y)) || throw(DimensionMismatch())
    return ProcessTensor(x.parent * y.parent)
end

# CanonicalMPO 无原生 dot：经 vectorize 转为 CanonicalMPS 视图求内积
# （per-site scaling 折入的约定两边一致）
LinearAlgebra.dot(ρA::ProcessTensor, ρB::ProcessTensor) =
    dot(vectorize(ρA.parent), vectorize(ρB.parent))



"""
    addition of two MPOs
"""
# 块对角直和：FiniteMPSAlgorithms 的 exact 算术限定 plain same-kind 链（MPO + MPO），
# 这里恢复 CanonicalMPO 的加法语义——两边 scaling 逐站点折入数据（scaling^L 约定），
# 直和后以 scaling = 1 重新包装；`-` 走 `x + (-y)` 自动恢复
function Base.:+(x::ProcessTensor, y::ProcessTensor)
    (length(x) == length(y)) || throw(DimensionMismatch())
    a = MPO(scaling(x.parent) .* x.parent.data)
    b = MPO(scaling(y.parent) .* y.parent.data)
    return ProcessTensor(CanonicalMPO(_plus_data(a, b).data))
end
# adding mpo with adjoint mpo will return an normal mpo
Base.:-(x::ProcessTensor, y::ProcessTensor) = x + (-y)
Base.:-(x::ProcessTensor) = -1 * x


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

# swap gate：委托给 FiniteMPSAlgorithms 的 CanonicalMPO `swap!`（Hastings
# 更新，与 TEMPO/GTEMPO 对齐；右规范形式与记录的键谱在截断误差内保持），
# 需要时自动完成规范化。`permute!` / `permute` 为链级置换的对应委托。
swap!(x::ProcessTensor, bond::Int; trunc::TruncationScheme=DefaultITruncation) = (swap!(x.parent, bond; trunc); x)

permute!(x::ProcessTensor, perm::AbstractVector{Int}; kwargs...) = (permute!(x.parent, perm; kwargs...); x)
permute(x::ProcessTensor, perm::AbstractVector{Int}; kwargs...) = ProcessTensor(permute(x.parent, perm; kwargs...))
