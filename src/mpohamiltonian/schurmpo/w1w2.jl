# arXiv:1407.1832v1 "Time-evolving a matrix product state with long-ranged interactions"
abstract type TimeEvoMPOAlgorithm <: MPSAlgorithm end
"""
    FirstOrderStepper

Abstract type for first-order time evolution stepping algorithms; concrete implementations include `WI` and `WII`.
"""
abstract type FirstOrderStepper <: TimeEvoMPOAlgorithm end
abstract type SecondOrderStepper <: TimeEvoMPOAlgorithm end

"""
    WI(tol, maxiter)
    WI(; tol=Defaults.tol, maxiter=Defaults.maxiter)

Parameters for the first-order W-type (WI) time evolution stepping algorithm: `tol` is the truncation error and `maxiter` the maximum number of iterations.
"""
struct WI <: FirstOrderStepper
	tol::Float64
	maxiter::Int
end
"""
    WI(; tol=Defaults.tol, maxiter=Defaults.maxiter)

Keyword constructor for `WI`; the default values are taken from `Defaults`.
"""
WI(; tol::Real=Defaults.tol, maxiter::Int=Defaults.maxiter) = WI(convert(Float64, tol), maxiter)

"""
    WII(tol, maxiter)
    WII(; tol=Defaults.tol, maxiter=Defaults.maxiter)

Parameters for the first-order W-type (WII) time evolution stepping algorithm: `tol` is the truncation error and `maxiter` the maximum number of iterations.
"""
struct WII <: FirstOrderStepper
	tol::Float64 
	maxiter::Int 
end
"""
    WII(; tol=Defaults.tol, maxiter=Defaults.maxiter)

Keyword constructor for `WII`; the default values are taken from `Defaults`.
"""
WII(; tol::Real=Defaults.tol, maxiter::Int=Defaults.maxiter) = WII(convert(Float64, tol), maxiter)

"""
    ComplexStepper(stepper)

Combine a first-order stepper `stepper` into a second-order complex-time stepper:
evolving twice with `dt₁=(1-im)*dt/2` and `dt₂=(1+im)*dt/2` yields second-order accuracy.
"""
struct ComplexStepper{F<:FirstOrderStepper} <: SecondOrderStepper
	stepper::F
end

function get_A(x::SchurMPOTensor)
	s1, s2 = size(x)
	# return x.Os[2:s1-1, 2:s2-1]
	return [x[i, j] for i in 2:s1-1, j in 2:s2-1]

end
function get_B(x::SchurMPOTensor)
	s1, s2 = size(x)
	# return x.Os[2:s1-1, end]
	return [x[i, s2] for i in 2:s1-1]
end
function get_C(x::SchurMPOTensor)
	s1, s2 = size(x)
	# return x.Os[1, 2:s2-1]
	return [x[1, j] for j in 2:s2-1]
end
get_D(x::SchurMPOTensor) = x[1, end]

function _SiteW_impl(WA, WB, WC, WD)
	s1, s2 = size(WA)
	r = Matrix{Any}(undef, s1+1, s2+1)
	r[1, 1] = WD
	for l in 2:s2+1
		r[1, l] = WC[l-1]
	end
	for l in 2:s1+1
		r[l, 1] = WB[l-1]
	end
	r[2:end, 2:end] .= WA
	return SparseMPOTensor(r)
end

function _sqrt2(dt::Complex) 
	r = sqrt(dt)
	return r, r
end

function _sqrt2(dt::Real)
	if dt >= zero(dt)
	 	r = sqrt(dt)
	 	return r, r
	 else
	 	r = sqrt(-dt)
	 	return r, -r
	end 
end

"""
    timeevompo(m, dt, alg)

Evolve a `SchurMPOTensor` or `MPOHamiltonian` by one step `dt` according to the time evolution stepping algorithm `alg`,
returning a new MPO tensor (or `MPOHamiltonian`). `alg` can be `WI`, `WII`, `ComplexStepper`,
or any `FirstOrderStepper`.
"""
function timeevompo(m::SchurMPOTensor, dt::Number, alg::WI)
	WA = get_A(m)
	δ₁, δ₂ = _sqrt2(dt)
	WB = get_B(m) .* δ₁
	WC = get_C(m) .* δ₂
	D = get_D(m)
	WD = _eye(scalartype(D), size(D, 1)) + dt * D
	return _SiteW_impl(WA, WB, WC, WD)
end

"""
    timeevompo(m::SchurMPOTensor, dt, alg::WII)

The `WII` stepper: constructs an MPO block containing the matrix exponential for each site.
"""
function timeevompo(m::SchurMPOTensor, dt::Number, alg::WII)
	A, B, C, D = get_A(m), get_B(m), get_C(m), get_D(m)
	Ddt = dt * D
	WD = exp(Ddt)
	d = phydim(m)
	mo = zero(Ddt)
	s1, s2 = size(A)
	δ₁, δ₂ = _sqrt2(dt)
	
	T = typeof(Ddt)

	WA = Array{T, 2}(undef, size(A))
	WB = Array{T, 1}(undef, size(B))
	WC = Array{T, 1}(undef, size(C))


	for a in 1:s1, b in 1:s2
		# println(size(Ddt), " ", size(mo), " ", size(C[b]), " ", size(B[a]), " ", size(A[a, b]))
		tmp = [Ddt mo mo mo; δ₂*C[b] Ddt mo mo; δ₁*B[a] mo Ddt mo; A[a,b] δ₁*B[a] δ₂*C[b] Ddt]
        m = exp(tmp)
        m = m[:, 1:d]
        WC[b] = m[(d+1):2*d, :]
        WB[a] = m[(2*d+1):3*d, :]
        WA[a, b] = m[(3*d+1):4*d, :]
	end
	return _SiteW_impl(WA, WB, WC, WD)
end

"""
    timeevompo(m::MPOHamiltonian{<:SchurMPOTensor}, dt, alg::FirstOrderStepper)

Evolve each site tensor of the `MPOHamiltonian` site by site.
"""
timeevompo(m::MPOHamiltonian{<:SchurMPOTensor}, dt::Number, alg::FirstOrderStepper) = MPOHamiltonian([timeevompo(mj, dt, alg) for mj in m.data])
"""
    timeevompo(h, dt, alg::ComplexStepper)

Evolve using the second-order complex-time stepper, returning the two evolution results `(U₁, U₂)`.
"""
function timeevompo(h::Union{SchurMPOTensor, MPOHamiltonian{<:SchurMPOTensor}}, dt::Number, alg::ComplexStepper)
	dt1, dt2 = complex_stepper(dt)
	return timeevompo(h, dt1, alg.stepper), timeevompo(h, dt2, alg.stepper)
end
"""
    timeevompo(m::MPOHamiltonian{<:SchurMPOTensor}, dt; alg=WII())

Keyword version of `timeevompo`, which defaults to `WII()`.
"""
timeevompo(m::MPOHamiltonian{<:SchurMPOTensor}, dt::Number; alg::TimeEvoMPOAlgorithm = WII()) = timeevompo(m, dt, alg)

"""
	complex_stepper(dt::Number)
Retun dt₁, dt₂
If U = exp(H*dt) is a first order stepper,
then U₁ = exp(H*dt₁), U₂ = exp(H*dt₂), and
U = U₁U₂ will be a second order stepper
"""
complex_stepper(dt::Number) = (1-im) * dt/2, (1+im) * dt/2
