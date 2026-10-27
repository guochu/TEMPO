using Test, LinearAlgebra, TensorOperations, ImpurityModelBase, QuAPI
using Random

using TEMPO
using TEMPO: QR, QRpos, LQ, LQpos, SVD, SDD, Polar, TimeAscending, TimeDscending, Zvalue2
import TEMPO: branches

using QuAPI.Dense
using QuAPI.Dense: BosonicSplitDelta
using QuAPI: BosonicContourIndex
using FiniteMPSAlgorithms: NoTruncation

Random.seed!(20260916)

include("util.jl")

# dense coefficient tensor of a ProcessTensor (MPO):
# axes = (ket_1, bra_1, ket_2, bra_2, ...), lattice site order
function _tempo_pt_array(x::TEMPO.ProcessTensor)
    L = length(x)
    d = size(x[1], 2)
    M = TEMPO.todense(x.parent; order=:lsb)
    T = reshape(M, ntuple(_ -> d, L)..., ntuple(_ -> d, L)...)
    return permutedims(T, vcat([[p, L + p] for p in 1:L]...))
end

# 只跑新增的 hyb × algorithm × contour 矩阵测试集（迭代用；正式入口在 api/quapi_crosscheck.jl）
source = read(joinpath(@__DIR__, "api", "quapi_crosscheck.jl"), String)
marker = "@testset \"QuAPI dense cross-check: hyb × algorithm × contour matrix\" begin"
idx = findfirst(marker, source)
@testset "newmatrix-only" begin
    include_string(Main, source[idx.start:end], "api/quapi_crosscheck.jl")
end
