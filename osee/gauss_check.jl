using TEMPO, ImpurityModelBase, LinearAlgebra

omic_spectrum(w, α, wc) = 2α * w * exp(-w / wc)
N, α, wc = 4, 0.1, 2.0
δτ = 0.1; β = N * δτ
spec = spectrum(w -> omic_spectrum(w, α, wc), lb=0, ub=wc)
bath = bosonicbath(spec, β=β)
lattice = PTLattice(N=N, δτ=δτ, d=3, contour=:imag)
Δ = Matrix(correlationfunction(bath, lattice).data)

function op_entropy(F, d, N, l)
    R = reshape(F, ntuple(_ -> d, 2N)...)
    rowdims = [N - s + 1 for s in 1:N]
    coldims = [N + (N - s + 1) for s in 1:N]
    perm = vcat(vcat(rowdims[1:l], coldims[1:l]), vcat(rowdims[l+1:N], coldims[l+1:N]))
    R2 = permutedims(R, perm)
    Mmat = reshape(R2, d^(2l), d^(2(N - l)))
    sv = svdvals(Mmat); p = sv .^ 2; p ./= sum(p); p = p[p .> 1e-18]
    return -sum(p .* log.(p))
end

function exact_osee(Δ, d; sgn=+1)
    N = size(Δ, 1)
    terms = [adaga(i, j, coeff=Δ[i, j]) for i in 1:N for j in 1:N if Δ[i, j] != 0]
    H = bosonoperator(quadratichamiltonian(N, terms); d=d)
    F = exp(sgn * H)
    return [op_entropy(F, d, N, l) for l in 1:N-1]
end

# ---------- 高斯协方差公式（正确的 c=(a,b†) 构造） ----------
function cov_sxsp(M)
    n = size(M, 1)
    W = zeros(2n, 2n); W[1:n, n+1:2n] .= M; W[n+1:2n, 1:n] .= transpose(M)
    G = inv(I - W * W)          # <c c†>
    K = G * W                   # <c c>
    sx = real.(G) + real.(K) - 0.5 * I
    sp = real.(G) - real.(K) - 0.5 * I
    return sx, sp
end

function entropy_sxsp(sx, sp, idx)
    m = length(idx)
    σ = zeros(2m, 2m)
    σ[1:m, 1:m] .= sx[idx, idx]
    σ[m+1:2m, m+1:2m] .= sp[idx, idx]
    Ω = zeros(2m, 2m); Ω[1:m, m+1:2m] .= Matrix(I, m, m); Ω[m+1:2m, 1:m] .= -Matrix(I, m, m)
    ν = real.(eigvals(im * Ω * σ))
    return sum(((v + 0.5) * log(v + 0.5) - (v - 0.5) * log(v - 0.5) for v in ν if v > 0.5 + 1e-12); init=0.0)
end

function gauss_osee(M, L)
    sx, sp = cov_sxsp(M)
    return entropy_sxsp(sx, sp, vcat(L, size(M, 1) .+ L))
end

# 单对压缩态自检：|ψ> = Σ s^n |n,n>, 只看一个 mode 的熵
analytic(s) = -log(1 - s^2) - (s^2 / (1 - s^2)) * log(s^2)
println("--- 单对压缩态自检 ---")
for s in (0.2, 0.5, 0.8)
    sx, sp = cov_sxsp(reshape([s], 1, 1))
    println("s=$s  gauss=", round(entropy_sxsp(sx, sp, [1]), digits=10),
            "  analytic=", round(analytic(s), digits=10))
end
# 两模真空（M=0）整体系熵应为 0
sx0, sp0 = cov_sxsp(zeros(2, 2))
println("vacuum whole=", round(entropy_sxsp(sx0, sp0, [1, 2]), digits=12))

# ---------- exp(-H): 有界 Gaussian 算符，OSEE 应收敛 ----------
println("\n--- exp(-H) 精确 quadratic，d 扫描 ---")
for d in (3, 4, 6)
    println("d=$d  S = ", round.(exact_osee(Δ, d; sgn=-1), digits=10))
end

println("\n--- 高斯协方差公式作用于 M = exp(-Δ) ---")
M = exp(-Δ)
println("sv(M) = ", round.(svdvals(M), digits=10))
for l in 1:N-1
    println("l=$l  S_gauss = ", round(gauss_osee(M, collect(1:l)), digits=10))
end

println("\n--- 对照：同一 M 的 off-diagonal block SVD 快速公式 ---")
for l in 1:N-1
    s = svdvals(Δ[1:l, l+1:N])
    println("l=$l  S_fast = ", round(2 * sum(x^2 * (1 - log(x^2)) for x in s if x > 0), digits=10))
end
