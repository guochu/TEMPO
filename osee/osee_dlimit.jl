# d -> ∞ 行为检验：
#   1) N=2，按总粒子数分块精确计算 exp(H) 的 OSEE 随 Fock 截断 d 的变化；
#   2) 高斯协方差公式给出 d -> ∞ 的极限值，与分块结果对照。
#
# 运行：JULIA_LIBSTDCXX_PROBE=0 julia --startup-file=no osee_dlimit.jl

using TEMPO, ImpurityModelBase, LinearAlgebra

omic_spectrum(w, α, wc) = 2α * w * exp(-w / wc)

N, α, wc = 2, 0.1, 2.0
δτ = 0.1
β = N * δτ

spec = spectrum(w -> omic_spectrum(w, α, wc), lb=0, ub=wc)
bath = bosonicbath(spec, β=β)
Δ = Matrix(correlationfunction(bath, PTLattice(N=N, δτ=δτ, d=2, contour=:imag)).data)
println("Δ = ", Δ)
println("eig(Δ) = ", eigvals(Δ))

# ---------- 1) 总粒子数分块精确计算 exp(H) 的 OSEE ----------
# H = Σᵢⱼ Δᵢⱼ a†ᵢaⱼ 守恒总粒子数 M = n₁+n₂，逐 M 块构造并取 SVD。
function block_entropy(Δ, d)
	ss = Float64[]
	for M in 0:2(d-1)
		n1s = collect(max(0, M - d + 1):min(d - 1, M))
		L = length(n1s)
		idx = Dict((n1s[k], M - n1s[k]) => k for k in 1:L)
		H = zeros(L, L)
		for k in 1:L
			n1, n2 = n1s[k], M - n1s[k]
			H[k, k] = Δ[1, 1] * n1 + Δ[2, 2] * n2
			if haskey(idx, (n1 - 1, n2 + 1))
				H[idx[(n1 - 1, n2 + 1)], k] += Δ[1, 2] * sqrt(n1 * (n2 + 1))
			end
			if haskey(idx, (n1 + 1, n2 - 1))
				H[idx[(n1 + 1, n2 - 1)], k] += Δ[2, 1] * sqrt((n1 + 1) * n2)
			end
		end
		append!(ss, svdvals(exp(H)))
	end
	p = ss .^ 2
	p ./= sum(p)
	p = p[p .> 1e-18]
	return -sum(p .* log.(p))
end

println("\n--- 分块精确 OSEE (N=2, l=1) ---")
for d in (4, 8, 16, 32, 64, 128, 192, 256)
	t = @elapsed S = block_entropy(Δ, d)
	println("d=$d  S=", round(S, digits=8), "  t=", round(t, digits=2), "s")
end

# ---------- 2) 高斯协方差公式：d -> ∞ 极限 ----------
# 取 M = exp(-Δ)；exp(+Δ) 与 exp(-Δ) 的归一化算符只差每个模的
# 粒子-空穴反射（局域操作），时间 cut 的纠缠熵相同。
function cov_sxsp(M)
	n = size(M, 1)
	W = zeros(2n, 2n)
	W[1:n, n+1:2n] .= M
	W[n+1:2n, 1:n] .= transpose(M)
	G = inv(I - W * W)              # <c c†>
	K = G * W                       # <c c>
	sx = real.(G) + real.(K) - 0.5 * I
	sp = real.(G) - real.(K) - 0.5 * I
	return sx, sp
end

function entropy_sxsp(sx, sp, idx)
	m = length(idx)
	σ = zeros(2m, 2m)
	σ[1:m, 1:m] .= sx[idx, idx]
	σ[m+1:2m, m+1:2m] .= sp[idx, idx]
	Ω = zeros(2m, 2m)
	Ω[1:m, m+1:2m] .= Matrix(I, m, m)
	Ω[m+1:2m, 1:m] .= -Matrix(I, m, m)
	ν = real.(eigvals(im * Ω * σ))
	return sum(((v + 0.5) * log(v + 0.5) - (v - 0.5) * log(v - 0.5)
				for v in ν if v > 0.5 + 1e-12); init=0.0)
end

println("\n--- 高斯协方差公式（d -> ∞ 极限，l=1） ---")
for (name, M) in (("exp(-Δ)", exp(-Δ)), ("exp(+Δ)", exp(Δ)))
	sx, sp = cov_sxsp(M)
	println("M=$name  S_inf = ", round(entropy_sxsp(sx, sp, [1, 3]), digits=8))
end
