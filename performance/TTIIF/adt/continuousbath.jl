# ============================================================================
# 性能对比基准（连续玻色子热库）：XTRGIF vs TDVPIF
#
# 模型：independent bosons（费米杂质 + 连续谱 Leggett 玻色浴）
#       浴谱 J(ω) = (α/2)(ω^d/ωc^(d-1)) e^{-ω/ωc}，取 d=3, ωc=1（默认 α=1）
#
# 参考解：ImpurityModelBase 解析解
#   - 虚时间 G(τ)：independentbosons_Gτ（相互作用热态，与虚时间轮廓一致）
#   - 实时间 G>(t)：independentbosons_greater(spec, t, ρ_0; ...)
#     （直积初态版本：杂质任意 ρ_0 ⊗ 浴热态，与纯实时格点的边界条件一致）
#
# 运行方式：
#   OMP_NUM_THREADS=1 julia --project=. performance/XTRGIF/continuousbath.jl
#
# 初态约定（与 TEMPO 格点定义一致）：
#   - 虚时间轮廓：轮廓自带 τ 支准备相互作用热态 → 解析解默认版本
#   - 实时间轮廓：初态为直积态 ρ_imp(最大混合) ⊗ 浴热态 → 解析解 ρ_0 版本
# ============================================================================

using TEMPO, ImpurityModelBase, LinearAlgebra, Printf

# ------------------------------ 模型参数 -----------------------------------
ϵ_d = 0.5                      # 杂质能级
nop = [1.0 0; 0 0]             # 杂质密度算符 n_d
hyb = AdditiveHyb([1.0, 0])    # 对角耦合：V = n_d

spec = Leggett(d=3, ωc=1)      # 连续谱玻色浴

chi = 50
trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
trunc2 = truncdimcutoff(D=2*chi, ϵ=1.0e-10, add_back=0)

# 杂质初态：最大混合态（对应 TEMPO 实时格点的默认边界条件 ρ₀ = ones(d)）
ρ_imp = Matrix(1.0 * I(2) / 2)

op1 = [0.0 0; 1 0]   # d
op2 = [0.0 1; 0 0]   # d†

# 两种待比较的影响泛函算法（Vector of Pairs 保持顺序）
algs = [
	"XTRGIF" => XTRGIF(k=5, fast=true, algmult=SVDCompression(trunc)),
	"TDVPIF"                 => TDVPIF(trunc=trunc, δ=0.1),
]

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
# 1. 虚时间基准：G(τ) = <d(τ) d†(0)>
# ============================================================================
println("="^70)
println("1. 虚时间基准 (imaginary time, Leggett 连续谱)")
println("="^70)

Nτ = 10
δτ = 0.1
β = Nτ * δτ

lattice = ADTLattice(N=Nτ, δτ=δτ, contour=:imag)
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

model = ImpurityHamiltonian(ϵ_d .* nop)

# 解析参考：相互作用热态 <d(τ) d†(0)>_β
exactGτ = independentbosons_Gτ(spec, β=β, ϵ_d=ϵ_d, Nτ=Nτ)

imag_results = []
mpsI_imag = Dict{String, Any}()

for (name, alg) in algs
	# --- 影响泛函构建 ---
	t_if = @elapsed mpsI = hybriddynamics(lattice, corr, hyb, alg)
	mpsI_imag[name] = mpsI
	Dmax = _maxbond(mpsI)

	# --- 观测量扫描（算符插入法） ---
	t_obs = @elapsed begin
		mpsK = sysdynamics(lattice, model, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice)
		adt = mult(mpsK, mpsI, trunc=trunc2)
		Zval = integrate(adt)

		corrs = ComplexF64[]
		c2 = ContourIndex(1)
		for i in 1:Nτ
			c1 = ContourIndex(i)
			ct = (i == 1) ? ContourOperator(c1, op1 * op2) : ContourOperator([c1, c2], [op1, op2])
			mpsK2 = sysdynamics(lattice, model, ct, trunc=trunc)
			mpsK2 = boundarycondition!(mpsK2, lattice)
			adt2 = mult!(mpsK2, mpsI, trunc=trunc)
			push!(corrs, integrate(adt2) / Zval)
		end
	end

	err = _relative_error(real.(corrs), exactGτ[1:length(corrs)])
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
# 2. 实时间基准：G>(t) = -i <d(t) d†(0)>
# ============================================================================
println()
println("="^70)
println("2. 实时间基准 (real time, Leggett 连续谱)")
println("="^70)

Nt = 10
δt = 0.02
tmax = Nt * δt

lattice = ADTLattice(N=Nt, δt=δt, contour=:real)
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

# 解析参考：直积初态（杂质 ρ_imp ⊗ 浴热态）下的 greater Green 函数，
# 对应 TEMPO 实时格点的默认边界条件（ρ₀ = 最大混合态）
exactG = [independentbosons_greater(spec, tj, ρ_imp, β=β, ϵ_d=ϵ_d) for tj in 0:δt:tmax]

real_results = []
mpsI_real = Dict{String, Any}()

for (name, alg) in algs
	# --- 影响泛函构建 ---
	t_if = @elapsed mpsI = hybriddynamics(lattice, corr, hyb, alg)
	mpsI_real[name] = mpsI
	Dmax = _maxbond(mpsI)

	# --- 观测量扫描（算符插入法，:+ 分支） ---
	t_obs = @elapsed begin
		mpsK = sysdynamics(lattice, model, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lattice)
		adt = mult(mpsK, mpsI, trunc=trunc2)
		Zval = integrate(adt)

		corrs = ComplexF64[]
		c2 = ContourIndex(1, branch=:+)
		for i in 1:Nt
			c1 = ContourIndex(i, branch=:+)
			ct = (i == 1) ? ContourOperator(c1, op1 * op2) : ContourOperator([c1, c2], [op1, op2])
			mpsK2 = sysdynamics(lattice, model, ct, trunc=trunc)
			mpsK2 = boundarycondition!(mpsK2, lattice)
			adt2 = mult!(mpsK2, mpsI, trunc=trunc)
			push!(corrs, integrate(adt2) / Zval)
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
