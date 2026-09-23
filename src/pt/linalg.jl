


# 精确（未压缩）算符链乘积：委托给 FiniteMPSAlgorithms 的 scaling-safe 乘法
# （CanonicalMPO × CanonicalMPO 的结果自带 scaling(x)·scaling(y)）
function Base.:*(x::ProcessTensor, y::ProcessTensor)
    @assert !isempty(x)
    (length(x) == length(y)) || throw(DimensionMismatch())
    return ProcessTensor(x.parent * y.parent)
end



"""
    addition of two MPOs
"""
# 块对角直和：FiniteMPSAlgorithms 的实现会把两边的 scaling 折入数据
Base.:+(x::ProcessTensor, y::ProcessTensor) = ProcessTensor(x.parent + y.parent)
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
