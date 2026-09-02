# ============================================================================
# 性能对比基准（PT 路径）：XTRGIF vs TDVPIF
#
# 模型：自由玻色子杂质模型（free bosons），一个玻色子杂质模式（频率 ϵ_d）
#       线性耦合到连续谱玻色子热库（Ohmic 型谱）：
#           H = ϵ_d a†a + Σₖ ωₖ bₖ†bₖ + Σₖ gₖ (a†bₖ + a bₖ†)
#       哈密顿量是二次型，Green 函数由 (N+1)×(N+1) 跳跃矩阵
#           A = [ϵ_d  gᵀ; g  diag(ω)]
#       的数值严格解给出（ImpurityModelBase.freebosons_Gτ / freebosons_greater_lesser）。
#
# 运行方式：
#   OMP_NUM_THREADS=1 julia --project=. performance/XTRGIF/pt/freebosons.jl
#
# 比较内容：
#   1. 效率：影响泛函构建时间、观测量扫描时间、最终键维
#   2. 精度：虚时间 Green 函数 G(τ) 与实时间 greater Green 函数 G>(t)
#      相对严格解的相对误差，以及两种算法所得 IF 的相互偏差
#
# 说明：
#   - 玻色子杂质的局域希尔伯特空间截断为 d（虚时间 d=10，实时间 d=4），
#     PTLattice 需显式传 d
#   - 耦合算符 a†（NonDiagonalHyb，保持二次型）
#   - 实时间初态为占据数 1 的 Fock 态 ρ_imp = |1⟩⟨1|，对应严格解的系数
#     关联矩阵 ρ₀ = diag(n_imp=1, n_B(ω₁), ..., n_B(ω_N))
# ============================================================================

using TEMPO, ImpurityModelBase, LinearAlgebra, Printf

# ------------------------------ 模型参数 -----------------------------------
omic_spectrum(w, α, wc) = 2α * w * exp(-w / wc)

ϵ_d = 1.0
α = 0.1
wc = 2.0

spec = spectrum(w -> omic_spectrum(w, α, wc), lb=0, ub=wc)

chi = 30
trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
trunc2 = truncdimcutoff(D=2*chi, ϵ=1.0e-10, add_back=0)

# 两种待比较的影响泛函算法（相同的 Prony 展开保证公平）
algexpan = OverDeterminedProny(n=20, tol=1.0e-8)
algs = [
	"XTRGIF" => XTRGIF(k=5, fast=true, algmult=DMRGMult1(trunc), algexpan=algexpan),
	"TDVPIF"                 => TDVPIF(trunc=trunc, δ=0.1, algexpan=algexpan),
]

# --------------------------- 严格解的辅助函数 -------------------------------
"""
	hoppingmatrix(bath, ϵ_d)

自由玻色子模型的 (N+1)×(N+1) 跳跃矩阵 A = [ϵ_d gᵀ; g diag(ω)]。
注意 `spectrumcouplings` 返回谱值 |gₖ|²，需开平方。
"""
function hoppingmatrix(bath, ϵ_d)
	ws, gs = frequencies(bath), spectrumcouplings(bath)
	n = length(ws)
	A = zeros(n + 1, n + 1)
	A[1, 1] = ϵ_d
	for k in 1:n
		A[1, 1+k] = gs[k]
		A[1+k, 1] = gs[k]
		A[1+k, 1+k] = ws[k]
	end
	return A
end

# ------------------------------ 公共函数 -----------------------------------
"""
	_relative_error(num, ref)

数值结果相对参考解的相对误差 norm(num - ref) / norm(ref)。
"""
_relative_error(num, ref) = norm(num - ref) / norm(ref)

"""
	_maxbond(mps)

MPS 的最大键维。
"""
_maxbond(mps) = maximum(bond_dimensions(mps))

# ============================================================================
# 1. 虚时间基准：G(τ) = <a(τ) a†(0)>（热平衡态）
# ============================================================================
println("="^70)
println("1. 虚时间基准 (imaginary time, PT, free bosons)")
println("="^70)

Nτ = 10
δτ = 0.1
β = Nτ * δτ
d = 10          # 杂质玻色子局域希尔伯特空间截断

lattice = PTLattice(N=Nτ, δτ=δτ, d=d, contour=:imag)
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

adag = bosonadagoperator(d=d)
nop = bosondensityoperator(d=d)
hyb = NonDiagonalHyb(adag)

model = ImpurityHamiltonian(ϵ_d .* nop)

# 观测量算符：op1·op2 = a·a† 对应 G(τ) = ⟨a(τ) a†(0)⟩
op1 = adag'
op2 = adag

# 严格解参考：连续谱离散化（δw=0.01）后对角化跳跃矩阵
b2 = discretebath(bath, δw=0.01)
A = hoppingmatrix(b2, ϵ_d)
exactGτ = freebosons_Gτ(A, collect(0:δτ:(Nτ-1)*δτ), 1, 1; β=β)

imag_results = []
mpsI_imag = Dict{String, Any}()

for (name, alg) in algs
	# --- 影响泛函构建 ---
	t_if = @elapsed mpsI = hybriddynamics(lattice, corr, hyb, alg)
	mpsI_imag[name] = mpsI
	Dmax = _maxbond(mpsI)

	# --- 观测量扫描（算符插入法） ---
	t_obs = @elapsed begin
		mps = mult!(sysdynamics(lattice, model, trunc=trunc), mpsI, trunc=trunc2)
		Zval = integrate(lattice, mps)

		corrs = ComplexF64[]
		c2 = ContourIndex(1)
		for i in 1:Nτ
			c1 = ContourIndex(i)
			ct = (i == 1) ? ContourOperator(c1, op1 * op2) : ContourOperator([c1, c2], [op1, op2])
			push!(corrs, integrate(lattice, apply!(ct, lattice, deepcopy(mps))) / Zval)
		end
	end

	err = _relative_error(corrs, exactGτ[1:length(corrs)])
	push!(imag_results, (name=name, t_if=t_if, t_obs=t_obs, Dmax=Dmax, err=err))
end

# 两种算法得到的 IF 互相比较
names = [n for (n, _) in algs]
dIF = distance(mpsI_imag[names[1]], mpsI_imag[names[2]]) / norm(mpsI_imag[names[2]])

@printf("%-24s %12s %14s %10s %14s\n", "算法", "IF构建(s)", "观测量扫描(s)", "最大键维", "G(τ)相对误差")
for r in imag_results
	@printf("%-24s %12.3f %14.3f %10d %14.3e\n", r.name, r.t_if, r.t_obs, r.Dmax, r.err)
end
@printf("\n两种算法 IF 的相对距离 (XTRG-IF vs TDVP): %.3e\n", dIF)

# ============================================================================
# 2. 实时间基准：G>(t) = -i <a(t) a†(0)>（直积初态，占据数 1）
# ============================================================================
println()
println("="^70)
println("2. 实时间基准 (real time, PT, free bosons)")
println("="^70)

Nt = 10
δt = 0.1
tmax = Nt * δt
β = 2.0
d = 4

lattice = PTLattice(N=Nt, δt=δt, d=d, contour=:real)
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

adag = bosonadagoperator(d=d)
nop = bosondensityoperator(d=d)
hyb = NonDiagonalHyb(adag)

model = ImpurityHamiltonian(ϵ_d .* nop)

op1 = adag'
op2 = adag

# 杂质初态：占据数 1 的 Fock 态，对应严格解的 ρ₀[1,1] = n_imp = 1
ρ_imp = bosonoccupationoperator(1, d=d)

# 严格解参考：系数关联矩阵 ρ₀ = diag(1, n_B(ω₁), ..., n_B(ω_N))
b2 = discretebath(bath, δw=0.01)
A = hoppingmatrix(b2, ϵ_d)
ws2 = frequencies(b2)
ρ₀ = Matrix(Diagonal([1.0; [boseeinstein(β, ωₖ) for ωₖ in ws2]]))
exactG = freebosons_greater_lesser(A, ρ₀, collect(0:δt:tmax), 1, 1)[1]

real_results = []
mpsI_real = Dict{String, Any}()

for (name, alg) in algs
	# --- 影响泛函构建 ---
	t_if = @elapsed mpsI = hybriddynamics(lattice, corr, hyb, alg)
	mpsI_real[name] = mpsI
	Dmax = _maxbond(mpsI)

	# --- 观测量扫描（算符插入法，:+ 分支） ---
	t_obs = @elapsed begin
		mps = mult!(sysdynamics(lattice, model, trunc=trunc), mpsI, trunc=trunc2)
		Zval = integrate(lattice, initialstate!(deepcopy(mps), lattice, ρ_imp))

		corrs = ComplexF64[]
		c2 = ContourIndex(1, branch=:+)
		for i in 1:Nt
			c1 = ContourIndex(i, branch=:+)
			ct = (i == 1) ? ContourOperator(c1, op1 * op2) : ContourOperator([c1, c2], [op1, op2])
			mps2 = apply!(ct, lattice, deepcopy(mps))
			push!(corrs, integrate(lattice, initialstate!(mps2, lattice, ρ_imp)) / Zval)
		end
	end

	corrs = -im .* corrs
	err = _relative_error(corrs, exactG[1:length(corrs)])
	push!(real_results, (name=name, t_if=t_if, t_obs=t_obs, Dmax=Dmax, err=err))
end

dIF = distance(mpsI_real[names[1]], mpsI_real[names[2]]) / norm(mpsI_real[names[2]])

@printf("%-24s %12s %14s %10s %14s\n", "算法", "IF构建(s)", "观测量扫描(s)", "最大键维", "G>(t)相对误差")
for r in real_results
	@printf("%-24s %12.3f %14.3f %10d %14.3e\n", r.name, r.t_if, r.t_obs, r.Dmax, r.err)
end
@printf("\n两种算法 IF 的相对距离 (XTRG-IF vs TDVP): %.3e\n", dIF)

println()
println("="^70)
println("基准完成")
println("="^70)
