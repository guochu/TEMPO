# Quickstart

This page presents complete workflows for four classes of typical problems; the code can be copied and run as-is. Which section to use depends on the **coupling form** and the **observable you want**:

| Problem | Coupling | Lattice | Example |
|---|---|---|---|
| [Standard spin-boson model (diagonal coupling, Rabi type, ADT)](@ref) | Diagonal | `ADTLattice` | `docs/tutorials/spinboson/rabitype.jl` |
| [JC-type off-diagonal coupling (PT framework)](@ref) | Off-diagonal (conjugate pair) | `PTLattice` | `docs/tutorials/spinboson/jctype.jl` |
| [Bosonic impurity: imaginary-time evolution and Matsubara correlations](@ref) | Diagonal or off-diagonal | `ADTLattice`/`PTLattice`, `contour=:imag` | `benchmark/independentbosons.jl` |
| [Time-dependent system-bath coupling](@ref) | Any | Any | `src/tdinfluencefunctional/` |

All computations follow the same template:

```text
1. trunc  = truncation scheme
2. lattice = ADTLattice/PTLattice (contour + number of steps + step size)
3. hyb    = AdditiveHyb / NonDiagonalHyb / ... (coupling operators)
4. spec/bath = spectrum + bosonicbath (bath spectrum and temperature)
5. corr   = correlationfunction(bath, lattice)
6. mpsI   = hybriddynamics(lattice, corr, hyb[, alg])    # influence functional
7. mpsK   = sysdynamics(lattice, model) + boundarycondition! (system dynamics + initial state)
8. cache  = environments(...); expectationvalue(...)      # observables
```

## Standard spin-boson model (diagonal coupling, Rabi type, ADT)

The system is a two-level spin with $\hat{A}=\hat{\sigma}_x/2$ (diagonal coupling), and the bath spectrum is sub-Ohmic.

```julia
using TEMPO, ImpurityModelBase, LinearAlgebra

# Truncation scheme: max bond dimension chi, singular value threshold ϵ
trunc = truncdimcutoff(D=chi, ϵ=1.0e-12, add_back=0)

# 1. Define the real-time contour lattice
lattice = ADTLattice(N=Nt, δt=δt, contour=:real)

# 2. Define system operators (σ_x for the system Hamiltonian; the diagonal entries of σ_z for the diagonal coupling operator A)
x = Matrix{ComplexF64}([0 1; 1 0])
z = Matrix{ComplexF64}([-1 0; 0 1])
zdiag = [z[i,i] for i in 1:size(z,1)]

# 3. Define the bath and hybridization style
hyb  = AdditiveHyb(zdiag)                       # diagonal coupling (note: pass the vector of diagonal entries)
spec = spectrum(w -> 2π*α * w^s * wc^(1-s), lb=0, ub=wc)   # sub-Ohmic spectrum
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

# 4. Build the influence functional (IF); the result is an ADT (MPS)
mpsI = hybriddynamics(lattice, corr, hyb, trunc=trunc)

# 5. Build the bare system dynamics and apply the boundary condition (initial state)
model = ImpurityHamiltonian(Δ .* x)             # system Hamiltonian H_S
mpsK  = sysdynamics(lattice, model, trunc=trunc)
mpsK  = boundarycondition!(mpsK, lattice, ρ₀=ρimp)

# 6. Precompute environments for observables
cache = environments(lattice, mpsK, mpsI)

# 7. Compute local observables (e.g., ⟨σ_z(t)⟩)
obs = ComplexF64[]
for i in 1:Nt
    pos = index(lattice, i, branch=:+)
    m   = ADTTerm(pos, zdiag)
    v   = expectationvalue(m, cache)            # already divided by Z, giving the normalized expectation value
    push!(obs, v)
end
```

!!! note "Diagonal coupling can also use XTRG-IF"
    The `hybriddynamics(lattice, corr, hyb, trunc=trunc)` call above takes the `PartialIF` path. For larger local dimensions (e.g., spin-1, d=3) or when more accurate compression is needed, the `XTRGIF` algorithm (see [Manual](@ref)) is usually faster and more memory-efficient; see the discussion in [Practice guide](@ref).

!!! tip "Measuring off-diagonal operators on the ADT path: operator insertion"
    Besides `ADTTerm`, on an ADT lattice off-diagonal operators and two-point correlation functions can also be measured by **inserting operators into the system dynamics** (see `test/models/rabimodel.jl` for an imaginary-time test):

    ```julia
    # Single-point operator (can be any matrix, including off-diagonal entries)
    ct = ContourOperator(ContourIndex(1), op)          # op is any matrix
    mpsK = sysdynamics(lattice, model, ct, trunc=trunc)

    # Two-point correlation ⟨op2(t_i) op1(t_1)⟩: insert one operator at each of the two positions
    ct = ContourOperator([ContourIndex(i), ContourIndex(1)], [op2, op1])
    mpsK = sysdynamics(lattice, model, ct, trunc=trunc)

    mpsK = boundarycondition!(mpsK, lattice, ρ₀=ρimp)  # real time requires ρ₀
    mps2 = mult!(mpsK, mpsI, trunc=trunc)
    v = integrate(mps2) / integrate(mps)               # already normalized
    ```

    This approach requires rebuilding `mpsK` and remultiplying the IF for each insertion, so it suits a small number of operators/correlation functions; for batch measurements of single-point diagonal quantities, `ADTTerm` + an `environments` cache is more efficient. On an ADT, diagonal two-point correlations can likewise be measured with the multi-point form `ADTTerm((pos2, pos1), (v2, v1))` (`apply!` + `integrate(mps2)/integrate(mps)`).

## JC-type off-diagonal coupling (PT framework)

The system couples to the bath through $\hat{A}=\hat{\sigma}_-/2$ (the conjugate pair $\hat{A}^\dagger,\hat{A}$). This is the **off-diagonal coupling** case of the paper, and it requires the PT framework together with the translationally invariant influence functional algorithm.

```julia
lattice = PTLattice(N=Nt, δt=δt, contour=:real)   # note: PTLattice

hyb  = NonDiagonalHyb(sp)                        # off-diagonal coupling: op*a + op'*a'
spec = spectrum(w -> 2π*α * w^s * wc^(1-s), lb=0, ub=wc)
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

# IF construction algorithm: translationally invariant IF (XTRG-style),
# using DMRG-type MPO-MPO multiplication + Prony exponential expansion
algmult  = DMRGMult1(trunc, maxiter=10)
algexpan = OverDeterminedProny(n=n, tol=1.0e-8, verbosity=2)
alg = XTRGIF(k=k, fast=true, algmult=algmult, algexpan=algexpan, verbosity=2)

mpsI = hybriddynamics(lattice, corr, hyb, alg)   # yields a ProcessTensor (MPO)

# System dynamics + initial state
model = ImpurityHamiltonian(Δ .* z)
mpsK  = sysdynamics(lattice, model, trunc=trunc)
mps   = mult!(mpsK, mpsI, trunc=trunc)

# Observables: PT real time requires passing the initial density matrix to the environments
cache = environments(lattice, mps, ρ₀=ρimp)

obs = ComplexF64[]
for i in 1:Nt
    ind = ContourIndex(i, :+)
    m   = ContourOperator(ind, x)               # any system operator x
    v   = expectationvalue(m, cache)
    push!(obs, v)
end
```

!!! warning "Two routes for measuring off-diagonal observables"
    `ADTTerm` on the ADT path only supports vectors of diagonal entries. To measure off-diagonal quantities such as $\langle\sigma_x\rangle$, you can: (i) switch to `PTLattice` + `NonDiagonalHyb` (physically equivalent to `AdditiveHyb` for real diagonal operators) and measure in batch with `ContourOperator` (recommended; see the example in the previous section); or (ii) stay on the ADT path and use **operator insertion** (see the tip above). See [Practice guide](@ref) for details.

## Bosonic impurity: imaginary-time evolution and Matsubara correlations

Here the impurity itself is a bosonic mode (with the local Hilbert space truncated to `d`), and the Matsubara Green's function $\langle \mathcal{T}_\tau \hat{a}(\tau)\hat{a}^\dagger(0)\rangle$ is computed on the imaginary-time contour.

```julia
lattice = ADTLattice(N=N, δτ=δτ, contour=:imag)      # imaginary-time contour
# or, for off-diagonal coupling: PTLattice(...)

a = bosonaoperator(d=d)                              # bosonic annihilation operator (ImpurityModelBase)
n = bosondensityoperator(d=d)

hyb  = AdditiveHyb(diag(n))                          # diagonal coupling
spec = Leggett(d=1, ωc=1)                            # or a custom spectrum(...)
bath = bosonicbath(spec, β=β)
corr = correlationfunction(bath, lattice)

mpsI  = hybriddynamics(lattice, corr, hyb, trunc=trunc)
model = ImpurityHamiltonian(ϵ_d .* n)                # H_S
mpsK  = sysdynamics(lattice, model, trunc=trunc)
Zval  = integrate(mpsK, mpsI)                        # partition function Z

# Two-point correlation function: insert operators into the system dynamics
op1, op2 = [0 0; 1 0], [0 1; 0 0]
c1, c2   = ContourIndex(1), ContourIndex(1)
ct       = ContourOperator([c1, c2], [op1, op2])
mpsK2    = sysdynamics(lattice, model, ct, trunc=trunc)
v        = integrate(mpsK2, mpsI) / Zval
```

For a bosonic impurity with **off-diagonal coupling** (e.g., `benchmark/bosonicimpurity.jl`):

```julia
lattice = PTLattice(N=N, δt=δt, d=d, contour=:real)
hyb  = NonDiagonalHyb(a')
alg  = XTRGIF(k=5, fast=true,
                             algmult=DMRGMult1(trunc, initguess=:rand),
                             algexpan=OverDeterminedProny(n=20, tol=1.0e-8))
mpsI = hybriddynamics(lattice, corr, hyb, alg)
```

## Time-dependent system-bath coupling

For time-dependent coupling (see `src/tdinfluencefunctional/`), use the time-dependent hybridization styles: `AdditiveTdHyb`, `NonAdditiveTdHyb`, `NonDiagonalTdHyb`. They accept a function `f(t)` describing the time dependence of the coupling strength:

```julia
hyb = NonDiagonalTdHyb(op, t -> f(t))   # or AdditiveTdHyb / NonAdditiveTdHyb
```
