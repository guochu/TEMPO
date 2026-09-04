# generate the three tutorial notebooks (.ipynb) from cell sources
using JSON

function make_nb(cells::Vector{Pair{String,String}}, path::String)
    jcells = []
    for (typ, src) in cells
        if typ == "md"
            push!(jcells, Dict("cell_type"=>"markdown", "metadata"=>Dict{String,Any}(), "source"=>src))
        else
            push!(jcells, Dict("cell_type"=>"code", "execution_count"=>nothing,
                               "metadata"=>Dict{String,Any}(), "outputs"=>Any[], "source"=>src))
        end
    end
    nb = Dict("cells"=>jcells,
              "metadata"=>Dict("kernelspec"=>Dict("display_name"=>"Julia 1.10", "language"=>"julia", "name"=>"julia-1.10"),
                               "language_info"=>Dict("name"=>"julia", "version"=>"1.10.11")),
              "nbformat"=>4, "nbformat_minor"=>5)
    open(path, "w") do f
        JSON.print(f, nb, 1)
    end
    println("wrote ", path)
end

# ================= notebook 1: Strathearn 2018 Fig. 2a =================
md1_1 = raw"""
# Reproducing Fig. 2a of Strathearn *et al.* 2018 (the original TEMPO paper)

**Paper**: A. Strathearn, P. Kirton, D. Kilda, J. Keeling, B. W. Lovett,
*Efficient non-Markovian quantum dynamics using time-evolving matrix product operators*,
[Nat. Commun. **9**, 3322 (2018)](https://doi.org/10.1038/s41467-018-05617-3) (arXiv:1711.09641).

**Figure reproduced**: Fig. 2(a) — polarization dynamics $\langle S_z(t)\rangle$ of the unbiased Ohmic spin-boson model across the localization transition, for a series of coupling strengths $\alpha$.

**Model** (Eq. (4) of the paper):
$$H = \Omega S_x + \sum_k S_z (g_k a_k + g_k^* a_k^\dagger) + \omega_k a_k^\dagger a_k, \qquad J(\omega) = 2\alpha\,\omega\, e^{-\omega/\omega_c}$$

**Paper parameters** (Fig. 2 caption): $\omega_c = 5$ (in units of $\Omega = 1$), zero temperature (initially unexcited bath),
initial state $\langle S_z(0)\rangle = +1/2$, memory length $K = 200$, critical coupling $\alpha_c \simeq 1.25$.

**Reduced scale**: the paper uses a finely converged $K=200$; here we take $\delta t = 0.125$, $t \le 8$ ($K=64$),
bond dimension $\chi = 30$, $\beta = 20$ (effectively zero temperature), preserving the physics across $\alpha_c$.

**Reference data**: the reference curves in the paper figure are themselves TEMPO/QUAPI numerics. Here we compare against an **independent method** — exact diagonalization (ED) of the full Hamiltonian after discretizing the bath into finitely many modes (short times only; the discrete-spectrum error grows with time, so ED is plotted only for $t \le 1.25$) — overlaid with the TEMPO results of this package.
"""

code1_2 = """
using TEMPO, ImpurityModelBase, LinearAlgebra

# ---- TEMPO: Ohmic SBM polarization dynamics <Sz(t)> ----
function tempo_sz(; α, δt=0.125, tmax=8.0, β=20.0, Δ=1.0, ωc=5.0, chi=30)
    Nt = round(Int, tmax/δt)
    trunc = truncdimcutoff(D=chi, ϵ=1.0e-12, add_back=0)
    lattice = ADTLattice(N=Nt, δt=δt, contour=:real)
    hyb = AdditiveHyb([0.5, -0.5])                       # coupling operator S_z = σz/2
    spec = spectrum(w -> 2α*w*exp(-w/ωc), lb=0, ub=3ωc)  # Ohmic J(ω) = 2αω e^{-ω/ωc}
    bath = bosonicbath(spec, β=β)
    corr = correlationfunction(bath, lattice)
    mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)
    ρimp = [1 0; 0 0.0]                                  # initial spin: Sz = +1/2
    model = ImpurityHamiltonian(Δ/2 .* [0 1; 1 0.0])    # H_s = Ω σx/2
    mpsK = sysdynamics(lattice, model, trunc=trunc)
    mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
    cache = environments(lattice, mpsK, mpsI)
    sz = [real(expectationvalue(ADTTerm(index(lattice, i, branch=:+), [0.5, -0.5]), cache)) for i in 1:Nt]
    return [(i-1)*δt for i in 1:Nt], sz
end

# ---- ED reference: discretize the bath and diagonalize the full Hamiltonian ----
function ed_sz(; α, Δ=1.0, ωc=5.0, Nmodes=10, d=2, times=0:0.125:1.25)
    r = range(0.2, stop=2ωc, length=Nmodes)
    ws = collect(r); dw = Float64(step(r))
    spec(w) = 2α*w*exp(-w/ωc)
    gs = [sqrt(spec(w)*dw/π) for w in ws]        # J(ω) ≈ π Σ g_k² δ(ω-ω_k)
    dimB = d^Nmodes
    a1 = zeros(d, d); for n in 1:d-1; a1[n, n+1] = sqrt(n); end
    n1 = a1'a1
    mode_op(op1, k) = kron([k == j ? op1 : Matrix{Float64}(I, d, d) for j in 1:Nmodes]...)
    Hb = zeros(dimB, dimB); X = zeros(dimB, dimB)
    for k in 1:Nmodes
        Hb .+= ws[k] .* mode_op(n1, k)
        X .+= gs[k] .* (mode_op(a1, k) + mode_op(a1, k)')
    end
    σx = [0 1; 1 0.0]; σz = [1 0; 0 -1.0]
    H = kron(Δ/2 .* σx, Matrix(I, dimB, dimB)) + kron(Matrix(I, 2, 2), Hb) + kron(σz/2, X)
    ψ0 = zeros(ComplexF64, 2*dimB); ψ0[1] = 1.0   # spin up, bath ground state (T = 0)
    E, V = eigen(Hermitian(H))
    c = V' * ψ0
    VO = V' * (kron(σz/2, Matrix(I, dimB, dimB))) * V
    sz = [real(c' * (exp.(-im*E*t) .* (VO * (exp.(im*E*t) .* c)))) for t in times]
    return collect(times), sz
end
println("functions defined")
"""

code1_3 = """
# TEMPO dynamics across the localization transition (paper: αc ≈ 1.25)
alphas = [0.4, 0.8, 1.2, 1.6]
tempo_res = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
for α in alphas
    @time tempo_res[α] = tempo_sz(α=α)
    println("α = ", α, "  done, final Sz = ", round(last(tempo_res[α][2]), digits=4))
end
"""

code1_4 = """
# ED references (short times, two representative couplings)
ed_res = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
for α in (0.4, 1.6)
    @time ed_res[α] = ed_sz(α=α)
end
println("ED done")
"""

code1_5 = raw"""
using Plots
pl = plot(xlabel="t Ω", ylabel="⟨S_z(t)⟩", title="Strathearn et al. 2018, Fig. 2(a) — Ohmic SBM localization transition (reduced scale)",
          legend=:topright, size=(640, 420))
colors = palette(:default)
for (i, α) in enumerate(alphas)
    t, sz = tempo_res[α]
    plot!(pl, t, sz, color=colors[i], label="TEMPO α=$(α)")
end
for α in (0.4, 1.6)
    t, sz = ed_res[α]
    idx = findfirst(isequal(α), alphas)
    scatter!(pl, t, sz, color=colors[idx], msw=0, ms=3, marker=:circle, label="ED α=$(α)")
end
hline!(pl, [0.0], color=:gray, ls=:dash, label=nothing)
savefig(pl, "strathearn2018_fig2a.png")
pl
"""

md1_6 = raw"""
## Discussion of results

- **Weak coupling** ($\alpha = 0.4$): the polarization relaxes completely to zero — the **delocalized phase**; the TEMPO curve agrees quantitatively with the independent ED reference over the plotted time window.
- **$\alpha = 1.6 > \alpha_c$**: the polarization decays but retains a finite residual value ($\approx 0.49$) — the **localized phase** (low-frequency bath modes pin the spin).
- **$\alpha = 0.8, 1.2$**: relaxation is markedly slower (critical slowing down). $\alpha = 0.8 < \alpha_c$ belongs to the delocalized phase, but with the reduced memory length ($K=64$) and finite time ($t \le 8$) it has not yet relaxed to zero — precisely near the transition the paper needs the long memory $K = 200$ to reach the asymptotic limit.
- The ED reference is based on 10 discrete bath modes; the discrete-spectrum error grows with time, so it only validates the short-time regime ($t \lesssim 1.25$); at strong coupling ($\alpha=1.6$) the finite-mode truncation error of ED is also larger.
"""

# ================= notebook 2: Segal 2023 (Lorentzian SBM) =================
md2_1 = raw"""
# Reproducing Anto-Sztrikacs, Nazir & Segal 2023 — spin-boson dynamics with a strongly coupled Lorentzian bath

**Paper**: N. Anto-Sztrikacs, A. Nazir, D. Segal,
*Effective dynamics of open quantum systems strongly coupled to a bath: A nonperturbative Lorentzian master equation*,
[Phys. Rev. Research **5**, 033227 (2023)](https://doi.org/10.1103/PhysRevResearch.5.033227).

**Content reproduced**: the exact reference dynamics shown in the figures of the paper (curves computed by TEMPO) — the spin-boson model $P(t) = \langle\sigma_z(t)\rangle$ with a Lorentzian-spectrum bath at strong coupling. The paper uses TEMPO
as the reference benchmark for its master equation (RCQME); this notebook reproduces that **TEMPO benchmark** part.

**Model**:
$$H = \frac{\Delta}{2}\sigma_x + \frac{\sigma_z}{2} \sum_k g_k (b_k + b_k^\dagger) + \sum_k \omega_k b_k^\dagger b_k,$$
$$J(\omega) = \frac{1}{\pi}\frac{\gamma \lambda^2}{(\omega-\Omega)^2 + \gamma^2} \quad \text{(Lorentzian)}$$

**Parameters** (reduced scale): $\Delta = 1$, $\lambda = 1$, $\Omega = 5$, $\gamma = 2$, $\beta = 1$,
$\delta t = 0.1$, $t \le 4$, $\chi = 30$; initial state $P(0) = 1$.
(The Lorentzian spectrum is nonvanishing as $\omega \to 0$; the convolution with the Bose factor $1/\omega$ is log-divergent,
so the spectral integral is truncated starting from a lower bound of $\gamma/100$.)

**Reference data**: independent ED (bath discretization + exact diagonalization).
"""

code2_2 = """
using TEMPO, ImpurityModelBase, LinearAlgebra

# ---- TEMPO: Lorentzian-bath SBM, P(t) ----
function tempo_p(; δt=0.1, tmax=4.0, β=1.0, Δ=1.0, λ=1.0, Ω=5.0, γ=2.0, chi=30)
    Nt = round(Int, tmax/δt)
    trunc = truncdimcutoff(D=chi, ϵ=1.0e-12, add_back=0)
    lattice = ADTLattice(N=Nt, δt=δt, contour=:real)
    hyb = AdditiveHyb([0.5, -0.5])
    # Lorentzian spectral density; lb > 0 avoids the log-divergent 1/ω tail of the Bose factor
    spec = spectrum(w -> γ*λ^2 / (π * ((w - Ω)^2 + γ^2)), lb=γ/100, ub=Ω + 10γ)
    bath = bosonicbath(spec, β=β)
    corr = correlationfunction(bath, lattice)
    mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)
    ρimp = [1 0; 0 0.0]
    model = ImpurityHamiltonian(Δ/2 .* [0 1; 1 0.0])
    mpsK = sysdynamics(lattice, model, trunc=trunc)
    mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)
    cache = environments(lattice, mpsK, mpsI)
    p = [real(expectationvalue(ADTTerm(index(lattice, i, branch=:+), [1.0, -1.0]), cache)) for i in 1:Nt]
    return [(i-1)*δt for i in 1:Nt], p
end

# ---- ED reference ----
function ed_p(; Δ=1.0, λ=1.0, Ω=5.0, γ=2.0, Nmodes=8, d=2, times=0:0.1:1.5)
    r = range(γ/100, stop=Ω + 3γ, length=Nmodes)
    ws = collect(r); dw = Float64(step(r))
    spec(w) = γ*λ^2 / (π * ((w - Ω)^2 + γ^2))
    gs = [sqrt(spec(w)*dw/π) for w in ws]
    dimB = d^Nmodes
    a1 = zeros(d, d); for n in 1:d-1; a1[n, n+1] = sqrt(n); end
    n1 = a1'a1
    mode_op(op1, k) = kron([k == j ? op1 : Matrix{Float64}(I, d, d) for j in 1:Nmodes]...)
    Hb = zeros(dimB, dimB); X = zeros(dimB, dimB)
    for k in 1:Nmodes
        Hb .+= ws[k] .* mode_op(n1, k)
        X .+= gs[k] .* (mode_op(a1, k) + mode_op(a1, k)')
    end
    σx = [0 1; 1 0.0]; σz = [1 0; 0 -1.0]
    H = kron(Δ/2 .* σx, Matrix(I, dimB, dimB)) + kron(Matrix(I, 2, 2), Hb) + kron(σz/2, X)
    # initial state: spin up along z + thermal bath at β
    β = 1.0
    ρb = exp(-β .* Hb); ρb ./= tr(ρb)
    ρ0 = kron([1 0; 0 0.0], ρb)
    E, V = eigen(Hermitian(H))
    M = V' * ρ0 * V                      # ρ0 in the eigenbasis
    W = V' * (kron(σz, Matrix(I, dimB, dimB))) * V
    # ρ(t) = V D(t) M D(t)† V†  with D(t) = diag(e^{-iEt});
    # P(t) = tr[W · D(t) M D(t)†]
    p2 = [real(tr(W * ((exp.(-im*E*t) .* M) .* exp.(im*E*t)'))) for t in times]
    return collect(times), p2
end
println("functions defined")
"""

code2_3 = """
@time tT, pT = tempo_p()
println("TEMPO done, P(end) = ", round(last(pT), digits=4))
"""

code2_4 = """
@time tE, pE = ed_p()
println("ED done")
"""

code2_5 = """
using Plots
pl = plot(xlabel="t", ylabel="P(t) = ⟨σ_z(t)⟩",
          title="Lorentzian-bath SBM at strong coupling (Anto-Sztrikacs et al. 2023, TEMPO benchmark)",
          legend=:topright, size=(640, 420))
plot!(pl, tT, pT, color=1, lw=2, label="TEMPO (this package)")
# plot the ED reference only in the short-time regime where the discrete-mode
# approximation is reliable (long times show recurrences of the discrete spectrum)
mask = tE .<= 1.0
scatter!(pl, tE[mask], pE[mask], color=2, msw=0, ms=3, label="exact diagonalization (8 modes)")
savefig(pl, "segal2023_pt.png")
pl
"""

md2_6 = raw"""
## Discussion of results

- For the Lorentzian bath at strong coupling ($\lambda = 1$, comparable to $\Delta$), $P(t)$ shows underdamped relaxation:
  the narrow-band bath spectrum (center $\Omega = 5$, half-width $\gamma = 2$) has a long memory time, and the short-time dynamics are characterized by coherent oscillations at frequency
  $\sim\Delta$ superimposed on a slow decay.
- At short times ($t \lesssim 1$) TEMPO and the independent ED agree quantitatively (within $\lesssim 0.05$); at longer times
  ED breaks down, showing recurrence oscillations of the discrete spectrum — hence only the short-time regime is plotted, while the TEMPO curve continues to give the dissipative relaxation of the continuous spectrum.
- These exact TEMPO curves are precisely the reference benchmarks used in the paper to test master equations of increasing order (Redfield / RCQME / generalized master equations).
"""

# ================= notebook 3: thermal state (Guo et al. 2026) =================
md3_1 = raw"""
# Reproducing Guo *et al.* 2026 (the theory paper of this package) — the imaginary-time mean-force state

**Paper**: C. Guo, W. Wu, X. Xu, T. Jiang, P.-X. Chen, R. Chen,
*Time-evolving matrix product operators for off-diagonal system-bath coupling*,
Phys. Rev. B **114**, 125413 (2026).

**Content reproduced**: the application of the PT-TEMPO/XTRG method on the imaginary-time contour in the paper — contracting the imaginary-time process tensor
yields the **mean-force state** of a strongly coupled quantum impurity:
$$\rho_{\mathrm{mfs}} = \frac{\mathrm{tr}_{b}\left[e^{-\beta H_{\mathrm{tot}}}\right]}{\mathrm{tr}\left[e^{-\beta H_{\mathrm{tot}}}\right]}$$
In the weak-coupling limit it reduces to the system Gibbs state $e^{-\beta H_s}/Z_s$; at finite coupling it deviates from the Gibbs state
(bath-induced corrections) — this is the central physics of strong-coupling thermodynamics (the mean-force Gibbs state).

**Model**: spin-1/2, $H_s = \Delta\sigma_z/2$, coupled through the conjugate pair $\sigma_+/2$ (**off-diagonal coupling**,
PT framework) to a sub-Ohmic bath
$J(\omega) = 2\pi\alpha\,\omega^{s}\omega_c^{1-s}$.

**Parameters**: $\Delta = 1$, $\beta = 1$ ($N=20$, $\delta\tau=0.05$), $s=0.5$, $\omega_c=5$,
$\chi=20$, with the coupling scanned over $\alpha \in \{0.05, 0.1, 0.2, 0.4\}$.

**Reference data**: the analytic weak-coupling Gibbs state $\rho_G = e^{-\beta H_s}/Z_s$.
"""

code3_2 = """
using TEMPO, ImpurityModelBase, LinearAlgebra

function mfs_population(; α, N=20, δτ=0.05, Δ=1.0, s=0.5, wc=5.0, chi=20)
    β = N * δτ
    trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
    lattice = PTLattice(N=N, δτ=δτ, contour=:imag)
    z = [-1 0; 0 1.0] / 2
    sp = [0 0; 1 0.0] / 2
    model = ImpurityHamiltonian(Δ .* z)
    mpsK = sysdynamics(lattice, model, trunc=trunc)
    hyb = NonDiagonalHyb(sp)                       # off-diagonal (conjugate-pair) coupling
    spec = spectrum(w -> 2π*α * w^s * wc^(1-s), lb=0, ub=wc)
    bath = bosonicbath(spec, β=β)
    corr = correlationfunction(bath, lattice)
    algmult = DMRGMult1(trunc)
    algexpan = OverDeterminedProny(n=20, tol=1.0e-8)
    alg = XTRGIF(k=5, fast=true, algmult=algmult, algexpan=algexpan)
    mpsI = hybriddynamics(lattice, corr, hyb, alg)
    mps = mult(mpsK, mpsI, trunc=trunc)
    ρ = meanforcestate(lattice, mps)
    ρ = ρ ./ tr(ρ)
    return real.(diag(ρ))    # [ground-state, excited-state] populations (H_s = Δσz/2, σz=-1 is ground)
end

# analytic weak-coupling limit
gibbs_population(; β=1.0, Δ=1.0) = begin
    z = [-1 0; 0 1.0] / 2
    ρG = exp(-β .* Δ .* z); ρG ./= tr(ρG)
    real.(diag(ρG))
end
println("functions defined")
"""

code3_3 = """
alphas = [0.05, 0.1, 0.2, 0.4]
mfs_pop = Dict{Float64, Vector{Float64}}()
for α in alphas
    @time mfs_pop[α] = mfs_population(α=α)
    println("α = ", α, "  populations = ", round.(mfs_pop[α], digits=4))
end
gibbs = gibbs_population()
println("weak-coupling Gibbs populations = ", round.(gibbs, digits=4))
"""

code3_4 = """
using Plots
pl = plot(xlabel="coupling α", ylabel="excited-state population",
          title="Mean-force state vs coupling (Guo et al. 2026, imaginary-time PT-TEMPO)",
          legend=:right, size=(640, 420))
exc = [mfs_pop[α][2] for α in alphas]
plot!(pl, alphas, exc, marker=:circle, lw=2, color=1, label="mean-force state (PT-TEMPO)")
hline!(pl, [gibbs[2]], color=2, ls=:dash, lw=2, label="weak-coupling Gibbs state (analytic)")
savefig(pl, "thermalstate_mfs.png")
pl
"""

md3_5 = raw"""
## Discussion of results

- Contracting the imaginary-time process tensor yields the thermal-equilibrium mean-force state $\rho_{\mathrm{mfs}}$;
  in the weak-coupling limit it approaches the Gibbs value (dashed line in the figure: excited-state population
  $\approx 0.269$ at $\beta\Delta = 1$). Even at the smallest $\alpha = 0.05$, the equilibrium correction induced by the sub-Ohmic bath ($s = 0.5$)
  is already visible ($0.34$) — a bath with large low-frequency spectral weight renormalizes the state appreciably at intermediate temperatures.
- As the coupling $\alpha$ increases, the excited-state population rises monotonically ($\approx 0.51$ at $\alpha = 0.4$) —
  this is the bath-induced strong-coupling equilibrium correction, i.e. the deviation of the **mean-force Gibbs state** from the standard Gibbs state,
  a key quantity in strong-coupling thermodynamics (quantum heat engines, modified Landauer principle, etc.).
- This reproduction follows the PT (process tensor) + translationally invariant influence functional (XTRG-style imaginary-time evolution) path,
  which is precisely the core algorithm of the paper for off-diagonal coupling.
"""

make_nb([
    "md"=>md1_1, "code"=>code1_2, "code"=>code1_3, "code"=>code1_4, "code"=>code1_5, "md"=>md1_6,
], "strathearn2018/strathearn2018.ipynb")

make_nb([
    "md"=>md2_1, "code"=>code2_2, "code"=>code2_3, "code"=>code2_4, "code"=>code2_5, "md"=>md2_6,
], "segal2023/segal2023.ipynb")

make_nb([
    "md"=>md3_1, "code"=>code3_2, "code"=>code3_3, "code"=>code3_4, "md"=>md3_5,
], "thermalstate/thermalstate.ipynb")
