# ============================================================================
# 性能对比基准（PT 路径）：XTRGIF vs TDVPIF
#
# 模型：自由玻色子杂质模型（independent bosons，费米杂质 + 单模玻色浴
#       DiracDelta），在 PT 格点上比较两种影响泛函算法。
#       参考：ImpurityModelBase 解析解
#       - 虚时间 G(τ)：independentbosons_Gτ（相互作用热态）
#       - 实时间 G>(t)：independentbosons_greater（直积初态 ρ_imp ⊗ 浴热态）
#
# 运行方式：
#   OMP_NUM_THREADS=1 julia --project=. performance/XTRGIF/pt/independentbosons.jl
#
# 比较内容：
#   1. 效率：影响泛函构建时间、观测量扫描时间、最终键维
#   2. 精度：虚时间 Green 函数 G(τ) 与实时间 greater Green 函数 G>(t)
#      相对解析解的相对误差，以及两种算法所得 IF 的相互偏差
#
# PT 路径说明：
#   - 耦合用 NonAdditiveHyb(n_d) 包装（PT 的 IF 算法签名要求 GeneralHybStyle）
#   - 虚时间无初态准备；实时间通过 initialstate! 设定杂质初态 ρ_imp，
#     对应解析解的直积初态约定（ρ_imp = 最大混合态）
# ============================================================================

using TEMPO, ImpurityModelBase, LinearAlgebra, Printf

# ------------------------------ 模型参数 -----------------------------------
ϵ_d = 0.5                      # 杂质能级
nop = [1.0 0; 0 0]             # 杂质密度算符 n_d
hyb = NonAdditiveHyb(nop)      # 对角耦合 V = n_d（PT 要求 GeneralHybStyle）

spec = DiracDelta(1)           # 单模玻色浴：ω₀ = 1，单位耦合

chi = 50
trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
trunc2 = truncdimcutoff(D=2*chi, ϵ=1.0e-10, add_back=0)

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
# 1. 虚时间基准：G(τ) = <d(τ) d†(0)>（相互作用热态）
# ============================================================================
println("="^70)
println("1. 虚时间基准 (imaginary time, PT)")
println("="^70)

Nτ = 10
δτ = 0.1
β = Nτ * δτ

lattice = PTLattice(N=Nτ, δτ=δτ, contour=:imag)
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

model = ImpurityHamiltonian(ϵ_d .* nop)

op1 = [0.0 0; 1 0]   # d
op2 = [0.0 1; 0 0]   # d†

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
# 2. 实时间基准：G>(t) = -i <d(t) d†(0)>（直积初态）
# ============================================================================
println()
println("="^70)
println("2. 实时间基准 (real time, PT)")
println("="^70)

Nt = 10
δt = 0.02
tmax = Nt * δt

lattice = PTLattice(N=Nt, δt=δt, contour=:real)
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

# 杂质初态：最大混合态，对应 TEMPO 实时格点的默认边界条件
ρ_imp = 0.5 .* one(nop)

# 解析参考：直积初态（ρ_imp ⊗ 浴热态）下的 greater Green 函数
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
