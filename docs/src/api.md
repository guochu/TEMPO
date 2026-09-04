# API reference

This page lists the documentation of `TEMPO`'s exported symbols grouped by functionality (generated automatically from the docstrings). Usage:

```julia
using TEMPO, ImpurityModelBase, LinearAlgebra
```

Notes:
- `branch` and `index` are provided by `QuAPI`, and `scalartype`, `space_l`, `space_r` by `TensorOperations`; all are reexported through TEMPO;
- The exponential-expansion functions (`OverDeterminedProny`, `DeterminedProny`, `exponential_expansion`, `expansion_error`) are provided by the `ExpExp` package and reexported; see the signature quick reference in the [Exponential expansion (`ExpExp`, reexported)](@ref) section at the end of this page;
- `phydim`, `spectrum`, `bosonicbath`, etc. are defined in `ImpurityModelBase`.

## Truncation schemes

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["tensorops/tensorops.jl", "tensorops/truncation.jl", "tensorops/distance.jl",
         "conversions.jl", "defaults.jl"]
```

## Tensor factorizations and matrix operations

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["tensorops/matrixalgebra.jl", "tensorops/tensorfactorizations.jl", "tensorops/distance.jl"]
```

## MPS/DMRG multiplication algorithms

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["algorithms.jl"]
```

## Contour indices

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["contourindices.jl"]
```

## MPO Hamiltonians

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["mpohamiltonian/def.jl", "mpohamiltonian/abstractmpotensor.jl",
         "mpohamiltonian/sparsempotensor.jl", "mpohamiltonian/schurmpotensor.jl",
         "mpohamiltonian/mpohamiltonian.jl",
         "mpohamiltonian/schurmpo/schurmpo.jl", "mpohamiltonian/schurmpo/longrange.jl",
         "mpohamiltonian/schurmpo/exponentialdecay.jl", "mpohamiltonian/schurmpo/generaldecay.jl",
         "mpohamiltonian/schurmpo/w1w2.jl"]
```

## Augmented density tensor (ADT)

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["adt/def.jl", "adt/adt.jl", "adt/abstractdefs.jl", "adt/orth.jl",
         "adt/mult/mult.jl", "adt/mult/svdmult.jl", "adt/mult/iterativemult.jl",
         "adt/integrate.jl", "adt/util.jl", "adt/linalg.jl"]
```

## Process tensors

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["pt/def.jl", "pt/orth.jl", "pt/mult/mult.jl", "pt/mult/svdmult.jl",
         "pt/mult/iterativemult.jl", "pt/util.jl", "pt/linalg.jl"]
```

## Local operator terms (ADTTerm / FockTerm)

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["adtterms.jl", "fockterms.jl", "contouroperators.jl"]
```

## Lattices

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["adtlattices/adtlattices.jl", "adtlattices/fockordering.jl",
         "adtlattices/realtime.jl", "adtlattices/imaginarytime.jl", "adtlattices/mixedtime.jl",
         "ptlattices/ptlattices.jl", "ptlattices/realtime.jl", "ptlattices/imaginarytime.jl",
         "ptlattices/mixedtime.jl", "ptlattices/integrate.jl"]
```

## Bath correlation functions

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["correlationfunction.jl"]
```

## Influence functionals (IF)

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["influencefunctional/influencefunctional.jl",
         "influencefunctional/partialif/partialif.jl", "influencefunctional/partialif/util.jl",
         "influencefunctional/partialif/realtime.jl", "influencefunctional/partialif/imaginarytime.jl",
         "influencefunctional/partialif/mixedtime.jl",
         "influencefunctional/ptpartialif/ptpartialif.jl", "influencefunctional/ptpartialif/util.jl",
         "influencefunctional/ptpartialif/realtime.jl", "influencefunctional/ptpartialif/imaginarytime.jl",
         "influencefunctional/ptpartialif/mixedtime.jl",
         "influencefunctional/ttiif/ttiif.jl",
         "influencefunctional/ttiif/adt/adt.jl", "influencefunctional/ttiif/adt/imag.jl",
         "influencefunctional/ttiif/adt/real.jl",
         "influencefunctional/ttiif/pt/pt.jl", "influencefunctional/ttiif/pt/imag.jl",
         "influencefunctional/ttiif/pt/real.jl",
         "influencefunctional/tdvpif/tdvpif.jl"]
```

## Time-dependent and BEC influence functionals

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["tdinfluencefunctional/tdinfluencefunctional.jl",
         "tdinfluencefunctional/partialif/partialif.jl", "tdinfluencefunctional/ttiif/ttiif.jl",
         "becinfluencefunctional/becinfluencefunctional.jl"]
```

## Boundary conditions and initial states

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["boundarycondition.jl"]
```

## System models

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["models/def.jl", "models/models.jl", "models/dissipative.jl",
         "models/unitary/unitary.jl", "models/unitary/adt.jl", "models/unitary/pt.jl"]
```

## Observables

```@autodocs; canonical = false
Modules = [TEMPO]
Pages = ["observables/adt/adt.jl", "observables/adt/envs.jl", "observables/adt/transfer.jl",
         "observables/adt/gf.jl",
         "observables/pt/pt.jl", "observables/pt/envs.jl", "observables/pt/transfer.jl",
         "observables/pt/mixedtransfer.jl",
         "observables/observables.jl",
         "observables/correlations.jl", "observables/heatcurrents.jl"]
```

## Exponential expansion (`ExpExp`, reexported)

The exponential expansion fits data $f(x)$ ($x = 1,\dots,N$) to a sum of exponentials $f(x) \approx \sum_i \alpha_i \beta_i^{x}$, used to represent hybridization functions in a form that is evolvable with MPOs.

```julia
OverDeterminedProny(; n::Int=10, tol::Real=1.0e-8, verbosity::Int=1, stepsize=1)
DeterminedProny(; n::Int=10, tol::Real=1.0e-8, verbosity::Int=1, stepsize=1)
#   n         maximum number of expansion terms (iteration cap)
#   tol       relative error threshold (tol*norm(f)); converges early once reached
#   verbosity ≥1 prints a non-convergence warning, ≥2 prints convergence info
#   stepsize  uniform sampling step size; chosen automatically when nothing

exponential_expansion(f::Vector{<:Number}, [L::Int,] alg=OverDeterminedProny())
#   → (αs, βs) such that f(x) ≈ Σᵢ αᵢ βᵢˣ; when a function is passed, k=1:L is sampled automatically

expansion_error(f, p) / expansion_error(f, coeffs, alphas)
#   2-norm error of the expansion at the sampled points
```

## Default values (`src/defaults.jl`)

| Constant | Value | Purpose |
|---|---|---|
| `DefaultTruncation` | `truncdimcutoff(D=100, ϵ=1e-14)` | General default |
| `DefaultITruncation` | `truncdimcutoff(D=200, ϵ=1e-10)` | Default for IF construction / `mult` |
| `DefaultKTruncation` | `truncdimcutoff(D=1000, ϵ=1e-10)` | Default for system dynamics |
| `DefaultIntegrationTruncation` | `truncdimcutoff(D=10000, ϵ=1e-12)` | Initial-state absorption |
| `DefaultMPOTruncation` | `truncdimcutoff(D=10000, ϵ=1e-12)` | MPO compression |
| `DefaultMultAlg` | `DMRGMult1(DefaultITruncation)` | Default algorithm for `mult` |
