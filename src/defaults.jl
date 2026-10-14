
#default settings
module Defaults

# 子模块不继承父模块（TEMPO）的绑定，所需的类型/函数需显式引入：
# 截断方案与 WII 来自 FiniteMPSAlgorithms，Prony 展开来自 ExpExp，
# DMRG1 是 TEMPO 本地定义（见 algorithms.jl）
using FiniteMPSAlgorithms: TruncationScheme, truncdimcutoff, truncrelerr, WII
using ExpExp: OverDeterminedProny
import ..DMRG1

	const maxiter = 100 # for DMRG iteration
	const D = 100 # default bond dimension 
	const tolgauge = 1e-14 # for MPS truncation
	const tol = 1e-12 # for DMRG iteration
	const tollanczos = 1.0e-10 # for lanczos eigensolver
	const tolexp = 1.0e-8 # for local eigen in DMRG
	const verbosity = 1

	alg_trunc_if(; D::Int=D, ϵ::Real=tolgauge, add_back::Int=0) = truncdimcutoff(D=D, ϵ=ϵ, add_back=add_back)
	alg_trunc_sys(; ϵ::Real=tolgauge) = truncrelerr(ϵ)
	alg_mult(; trunc::TruncationScheme=alg_trunc_if()) = DMRG1(trunc=trunc)
	alg_evo() = WII()
	alg_expan(; n::Int=15, tol::Real=1.0e-4, verbosity::Int=0) = OverDeterminedProny(n=n, tol=tol, verbosity=verbosity)
end
