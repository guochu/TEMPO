# Internals

This document describes the internals of the `TEMPO` toolkit in depth: data structures, algorithm workflows, tensor conventions and numerical details, mapped one-to-one onto the formulas of the paper (Phys. Rev. B 114, 125413 (2026)). For interfaces and usage, see [Quickstart](@ref) and [Manual](@ref).

## 1. Overview and module structure

The toolkit implements the TEMPO algorithm on top of `ImpurityModelBase` (system model definitions) and `QuAPI` (basic tensor-network types `Dense1DTN`, `DenseMPSTensor`, `ContourIndex`, `branch`/`index`, etc.). The loading order is defined in `src/TEMPO.jl`:

```
tensorops ──→ defaults ──→ mpohamiltonian ──→ contourindices
      └─→ adt ──→ pt ──→ conversions ──→ adtterms / fockterms
      └─→ adtlattices ──→ ptlattices ──→ contouroperators
      └─→ correlationfunction ──→ influencefunctional ──→ tdinfluencefunctional
      └─→ boundarycondition ──→ models ──→ observables
```

Responsibilities of the individual modules:

| Module | Responsibility |
|---|---|
| `tensorops/` + `algorithms.jl` | Truncation schemes, exponential expansion, MPS algorithms (`DMRGMult1`, `SVDCompression`), `transfermatrix`, tensor utilities |
| `defaults.jl` | All default hyperparameters (`Defaults`) |
| `mpohamiltonian/` | `SchurMPOTensor`, `MPOHamiltonian`, long-range decay terms, `timeevompo` (WI/WII/ComplexStepper) |
| `adt/` | ADT (MPS) types, orthogonalization, SVD compression, `mult!`, integration/swap gates |
| `pt/` | `ProcessTensor` (MPO) types, orthogonalization, multiplication |
| `adtlattices/` | Real/imaginary/mixed-time ADT lattices and index mappings, Fock orderings |
| `ptlattices/` | Real/imaginary/mixed-time PT lattices, `integrate` (PT evaluation) |
| `correlationfunction.jl` | Bath correlation function wrappers (`IndexCorrelationFunction`, branch correlations) |
| `influencefunctional/` | Influence functional construction: `PartialIF` (bond dimension 2), `XTRGIF` (XTRG style), `TDVPIF` (second-order single-site TDVP imaginary-time flow), `PTPartialIF`, XTRG-IF for PT |
| `tdinfluencefunctional/` | Influence functionals for time-dependent coupling (`*TdHyb`) |
| `boundarycondition.jl` | Initial state/boundary conditions `boundarycondition!`, `initialstate!` |
| `models/` | System Hamiltonians/Lindblad operators, `sysdynamics` |
| `observables/` | Environment caches, expectation values, partition function, transfer matrices |

## 2. Core data structures

### 2.1 ADT (Augmented Density Tensor / MPS)

Defined in `src/adt/def.jl`:

```julia
struct ADT{T<:Number, R<:Real} <: Dense1DTN{T}
    data::Vector{Array{T, 3}}            # third-order site tensor
    s::Vector{Union{Missing, Vector{R}}} # Schmidt vectors (singular vectors)
    scaling::Ref{Float64}                # global scaling factor
end
```

- Tensor convention (three legs): leg 1 = left auxiliary (bond) index, with dimension 1 at the leftmost site; leg 2 = physical index (1 corresponds to $|0\rangle$, 2 corresponds to $|1\rangle$); leg 3 = right auxiliary index, with dimension 1 at the rightmost site.
- `s` holds the vector of singular values for each bond (length = number of sites + 1), used for the `iscanonical` check and the `easy_swap!` swap gate.
- `scaling` records the overall scaling (see §2.3).

### 2.2 ProcessTensor (MPO)

A `ProcessTensor` is a column of four-leg tensors (`Vector{Array{T,4}}`) with leg order `(left bond, physical in, physical out, right bond)`, representing a chain of influence-functional operators (an MPO). An `ADT` can be viewed as a special case of an MPO (the out leg contracted into the in leg, i.e., a "diagonalized" state).

### 2.3 The scaling mechanism

For numerical stability, all orthogonalization/truncation steps avoid normalizing the tensors directly; instead the norms are accumulated into the global `scaling` (`src/adt/orth.jl`):

- `_renormalize!(psi, r, normalize)`: computes the norm `nr` of `r`; if `!normalize`, calls `_rescaling!` to multiply `nr^(1/L)` into `scaling`, and then applies `lmul!(1/nr, r)`;
- `_rescaling!(psi, n)`: `scaling *= n^(1/L)`.

This keeps the bond tensors at O(1) magnitude and avoids overflow/underflow from exponential growth/decay. The result of `mult!` carries the scaling `scaling(x) * scaling(y)`.

## 3. Contours and lattices

### 3.1 Time contours and branches

- Imaginary time (Matsubara): branch `:τ`, interval $[0, \beta]$, scalar type `Float64`;
- Real time (Keldysh): branches `:+` and `:-`, interval $[0, t]$, scalar type `ComplexF64`;
- Mixed (Kadanoff-Baym): the imaginary-time branch plus the upper and lower real-time branches.

`branch`/`index` are provided by `QuAPI` and imported (`src/TEMPO.jl`). Each lattice is also associated with a `TimeOrderingStyle` (`TimeAscending`/`TimeDscending`) and a `LayoutStyle` (currently only `TimeLocalLayout`, i.e., the state is locally arranged at each time step, which makes it easy to apply time-local operators).

### 3.2 Fock orderings (FockOrdering)

Fock orderings define the arrangement conventions of the Grassmann operators on a lattice (`src/adtlattices/fockordering.jl`):

| Type | Contour | Meaning |
|---|---|---|
| `M2M1` (alias `MM`) | Imaginary time | Each imaginary-time step contains the two orderings M2, M1; the imaginary branch is in descending time order |
| `M2m2M1m1` (alias `MmMm`) | Real time | Each real-time step contains the two orderings and their conjugates m2, m1; the real branch is in descending time order |
| `M2M1_m1M1m2M2` | Mixed | Real branch in ascending time order, imaginary branch in descending time order |

### 3.3 Lattice types and index mappings

Real-time ADT (`src/adtlattices/realtime.jl`):

```julia
struct RealADTLattice1Order{O<:RealFockOrdering} <: RealADTLattice{O}
    δt::Float64; d::Int; N::Int; ordering::O
end
```

Derived properties: `t = N*δt`, `Nt = N`, `k = kt = N+1`, `length = 2k`. Index mapping (descending, with `+1` reserving room for the boundary Grassmann numbers):

```julia
index(i, branch=:+) = TL - 2i + 1    # TL = 2(N+1)
index(i, branch=:-) = TL - 2i + 2
```

That is, the physical sites are arranged as `(N,+) (N,−) (N−1,+) (N−1,−) … (1,+) (1,−)`, with positions 1 and 2 reserved for the initial-state boundary.

Imaginary-time ADT: `length = k = N+1`, `index(i) = k+1-i` (descending), with the vacuum boundary at the rightmost end (i=1).

PT lattices (`src/ptlattices/`):
- `ImagPTLattice1Order`: `length = N`, `index(i) = N-i+1`, only the `:τ` branch, scalar type `Float64`;
- `RealPTLattice1Order`: `length = 2N`, `index(i, :+) = 2(N-i)+1`, `index(i, :-) = 2(N-i)+2` (descending, with `index(N, :+)=1` and `index(1, :-)=2N`);
- `MixedPTLattice1Order`: `Nt` real-time steps + `Nτ` imaginary-time steps, with `index` interleaved between the real and imaginary branches.

`indexmappings(lattice)` returns a Dict{(timestep, branch) => global_index}, and `vacuumstate(T, lattice)` constructs an all-ones vacuum ADT.

### 3.4 Index-shift convention for mixed contours (important)

On the mixed (Kadanoff-Baym) contour, the mapping from "correlation-function/operator-insertion indices → lattice indices" is **different** on the ADT and PT paths; this is the detail most prone to errors in the implementation:

**ADT path (`MixedADTLattice1Order`, `src/adtlattices/mixedtime.jl`) — shifted by one on all branches.** The lattice stores `kτ = Nτ + 1` imaginary-time points and `kt = Nt + 1` real-time points per branch: lattice index 1 of each branch is the boundary/junction point (`τ=0` glued to the start of the `:-` branch, `τ=β` glued to the start of the `:+` branch, `t=0` glued to the end of the imaginary-time branch), handled by `boundarycondition!` and **receiving no influence-functional gates at all**; the correlation-function index `i` is mapped to the lattice index `i+1` uniformly on **all** branches (`:τ`, `:+`, `:-`), i.e., the lattice indices `2..kτ` (`2..kt`) carry the correlation indices `1..Nτ` (`1..Nt`). This convention is implemented consistently in `partialif`, `partialif_naive` (`src/influencefunctional/partialif/util.jl`) and `hybriddynamics!` (`src/influencefunctional/partialif/mixedtime.jl`) (in the form `i′ = (b == :τ || lattice isa MixedADTLattice) ? i+1 : i`). For comparison with the pure contours: the imaginary-time ADT lattice shifts only the `:τ` branch by one; the real-time ADT lattice maps directly without shifting.

**PT path (`MixedPTLattice1Order`, `src/ptlattices/mixedtime.jl`) — direct mapping, no shifting.** Unlike the mixed ADT lattice, the mixed PT lattice has no extra boundary points; it stores exactly `Nτ` imaginary-time points and `Nt` real-time points per branch (`length = 2Nt + Nτ`). Therefore `ContourIndex(i, branch)` directly addresses the lattice site at contour time `(i-1)δ` on the corresponding branch — the same direct mapping as for the pure imaginary/real-time PT lattices; `partialif_naive` (`src/influencefunctional/ptpartialif/`) and the `apply!` of `ContourOperator` are both implemented this way.

**Time alignment of operator insertions.** On either path, the `ContourIndex(i, branch)` of an inserted operator corresponds to the contour time `(i-1)δ` (`(i-1)δτ` on the `:τ` branch, `(i-1)δt` on the real branches): when the system dynamics applies the propagators row by row, `_get_U_prime` detects `ContourIndex(i, branch)` on the propagator of step `i` and multiplies the operator in. Therefore, when comparing Green's functions against an ED reference, the results for `ContourIndex(i)` must be aligned with the reference values at time `(i-1)δ`. For the verification of imaginary-time/real-time/mixed Green's functions with a single-mode bath + ED, see the mixed-time test sets in `test/models/rabimodel.jl` (ADT) and `test/ptmodels/toymodel.jl` (PT).

## 4. Bath correlation functions and exponential expansion

`correlationfunction.jl` wraps correlation functions into `IndexCorrelationFunction` (containing the ηᵢⱼ matrix of a `CorrelationMatrix`, generated by `exponential_expansion`); in real time, `branch(corr, :+, :-)` and the like extract the four sign-labeled components η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻.

Exponential expansion (the `ExpExp` package, re-exported by TEMPO):
- `OverDeterminedProny`: expands the correlation function into a finite sum of exponentials (least-squares Prony method);
- `DeterminedProny`: a deterministic variant of the Prony method;
- `expand_decayterm(decayterm, alg=...)` returns the coefficient list `(η₁, η₂, …)`;
- `expansion_error` estimates the expansion error.

## 5. Influence functionals for diagonal coupling (PartialIF / ADT path)

### 5.1 `partialif_densemps`: partial influence functional with bond dimension 2

The algorithm comes from Strathearn et al. (2018) and is implemented in `src/influencefunctional/partialif/util.jl`:

- Inputs: row position `row`, column positions `cols`, diagonal operator `op` (`z`), coefficients `coefs` (η values);
- Builds a single-row MPS (bond dimension = `d`, which is 2 for diagonal coupling with d=2):
  - at the row site, places the "one-body" tensor `exp(ηₛₛ z²)` (a diagonal matrix, `tmp[i,i,i]`);
  - at the remaining sites, places the "two-body" gates `exp(ηᵢⱼ z⊗z)` (`tmp[i,:,i] = m[i,:]`, with the bond index carrying the physical states);
  - the first/last sites are the boundary tensors `(1,d,d)`/`(d,d,1)`;
- `_fit_to_full` embeds the MPS into the whole lattice: the sites before and after the row are padded with all-ones vacuum tensors `ones(1,d,1)`, and the unreached positions within the row interval are filled with identity tensors of the `Σᵢ δ` type (`tmp[:,i,:] = I`).

### 5.2 `partialif_naive`

The naive version: for each column index in turn, `apply!` the exponential gate `exp(coef * op⊗op)` and then `canonicalize!` with truncation, finally obtaining the partial IF (same file). For `coef`, real time uses ηᵢⱼ (including its sign), while imaginary time uses ηᵢⱼ/2, split equally between the two halves.

### 5.3 `hybriddynamics!`

(`src/influencefunctional/partialif/realtime.jl` and the corresponding imaginary/mixed-time implementations) builds the partial IF row by row and multiplies them in one by one:

```julia
for i in 1:Nt, b1 in branches(lattice)
    pos1 = index(lattice, i, branch=b1)
    # collect the coefficients of this row against all columns
    tmp = partialif_densemps(ds, pos1, pos2s, op, coefs)
    mult!(gmps, tmp, trunc=trunc)
end
```

The partial IF of each `(i, b1)` row is multiplied into the global `gmps` and compressed by SVD (`DefaultITruncation`). This is the ADT implementation of Eqs. (2)-(8) of the paper. Note: the snippet above is the pure real-time lattice version (direct mapping); imaginary-time lattices (for the `:τ` branch) and mixed lattices (`MixedADTLattice1Order`, for **all** branches) must shift the correlation index `i` to the lattice index `i+1` (`pos1 = index(lattice, i+1, branch=b1)`); see §3.4 for details.

## 6. Translationally invariant influence functional (XTRG-IF)

This corresponds to the "translationally invariant + exponential expansion + MPO time evolution" scheme in the appendix of the paper; the entry points are `influenceoperator`/`influenceoperatorexponential`/`differentialinfluencefunctional` (`src/influencefunctional/ttiif/`).

### 6.1 Exponential expansion and `SchurMPOTensor`

`adt_ti_mpotensor` (`src/influencefunctional/ttiif/adt/imag.jl`):

```julia
m1 = GenericDecayTerm(op1, op2, corr.ηⱼₖ[2:end])   # long-range (cross-step) coupling
m2 = GenericDecayTerm(op2, op1, corr.ηₖⱼ[2:end])
m1s = expand_decayterm(m1, alg)                     # Prony expansion
m2s = expand_decayterm(m2, alg)
h1  = (corr.ηₖⱼ[1] + corr.ηⱼₖ[1]) .* (op1 * op2)    # same-time (diagonal block) term
return SchurMPOTensor(h1, vcat(m1s, m2s))
```

`SchurMPOTensor` is a block-triangular operator structure (`src/mpohamiltonian/schurmpo/`): `D = h1` (diagonal block, same-time coupling), `A = {decay terms}` (diagonal block, long-range coupling), and `B`, `C` are the upper/lower triangular connector blocks.

### 6.2 Time evolution: WI / WII / ComplexStepper

`timeevompo(m, dt, alg)` (`src/mpohamiltonian/schurmpo/w1w2.jl`); the schemes come from Zaletel et al., arXiv:1407.1832:

- **WI** (first order): `WD = I + dt·D`, `WB = B·√δt`, `WC = C·√δt` (`_sqrt2` returns `(√|dt|, −√|dt|)` for negative dt), which are then assembled into sparse MPO tensors;
- **WII** (first order, more accurate): construct the `4d×4d` block matrix

  ```
  [[Ddt,  0,   0,   0 ],
   [√δ₂C, Ddt, 0,   0 ],
   [√δ₁B, 0,   Ddt, 0 ],
   [A,    √δ₁B,√δ₂C,Ddt]]
  ```

  after taking `exp`, extract the blocks `WC`, `WB`, `WA` column by column, and `WD = exp(Ddt)`;
- **ComplexStepper**: `dt₁ = (1−i)dt/2`, `dt₂ = (1+i)dt/2`; two first-order steps are composed into a second-order accurate scheme, returning the two sets of MPOs before/after evolution (`timeevompo` returns `(U₁, U₂)`).

### 6.3 Fitting to the lattice: `_fit_to_lattice`

For the imaginary-time ADT, `_fit_to_lattice` (`src/influencefunctional/ttiif/adt/imag.jl`) tiles the 3 MPO tensors onto `L = N+1` sites in a translationally invariant pattern:

- position 1 (j=N+1) ← `mpstensors[1]` (left boundary tensor);
- positions 2…N−1 (j=N−1…3) ← `mpstensors[2]` (bulk tensor, repeated);
- position N (j=2) ← `mpstensors[3]` (right boundary tensor);
- position N+1 (j=1) ← `ones(1, d, 1)` (vacuum).

`_tompsj` contracts the out leg of the MPO tensor with an all-ones vector (`mpoj[1,2,3,4]*a[4]`), turning the operator into a "state" tensor so that it can be embedded into the ADT. The real-time ADT is similar but returns a tuple of the 4 branch MPOs.

### 6.4 Real-time branch MPOs and differential influence functionals

The PT real-time `influenceoperator` (`src/influencefunctional/ttiif/pt/real.jl`) returns the 4 branch MPOs `(η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻)`; `influenceoperatorexponential` first applies `timeevompo` to each branch and then fits it, with `FirstOrderStepper` returning 4 and `ComplexStepper` returning 8. `differentialinfluencefunctional` multiplies them successively (in the PT branch order `h2*h1`, `h3*…`, `h4*…`) to obtain the complete differential influence functional of a single time step.

## 7. The PT framework for off-diagonal coupling

This corresponds to the "effective Hamiltonian + copy tensor" scheme for off-diagonal coupling in the paper.

### 7.1 The effective Hamiltonian H_eff

Off-diagonal coupling between the system and the bath (e.g., the JC model) is brought to a diagonal form by introducing auxiliary "copy" degrees of freedom: in `H_eff`, the copied system states are paired with the bath operators, so that the influence functional remains a product of exponential gates. The toolkit represents the influence functional of this path with a PT (MPO).

### 7.2 `pt_ti_mpotensor` and symbolized correlation functions

An MPO tensor is constructed for each branch combination:

```julia
mpoj1 = pt_ti_mpotensor(η⁺⁺, op1, op2, :+, :+, algexpan)
mpoj2 = pt_ti_mpotensor(η⁺⁻, fused_op(op1, :+), fused_op(op2, :-), :+, :-, algexpan)
mpoj3 = pt_ti_mpotensor(η⁻⁺, fused_op(op1, :-), fused_op(op2, :+), :-, :+, algexpan)
mpoj4 = pt_ti_mpotensor(η⁻⁻, transpose(op1), transpose(op2), :-, :-, algexpan)
```

- the diagonal branches (++, −−) use the original operators directly (the −− branch takes the transpose, corresponding to the conjugate propagator);
- the off-diagonal branches (+−, −+) use `fused_op` to build "fused" two-leg operators: on the `:+` branch `op1` acts and the identity acts on the − branch; on the `:-` branch the identity acts and `transpose(op1)` acts on the − branch.

### 7.3 `fused_op` and `split_mpotensor`

- `fused_op(op, f)` (`src/influencefunctional/ttiif/pt/real.jl`): uses the identity-embedding tensor `f = reshape(I_{d²}, d², d, d)` to lift a one-leg operator to a two-leg `(d²×d²)` operator;
- `split_mpotensor(mpoj, trunc)`: reshapes the 4-leg MPO tensor into a 6-leg tensor with physical dimension `d²`, and splits it by SVD into the two factors `u` and `v` (both carrying `√s`), used to place the coupling of an off-diagonal branch at different positions on the two branches.

### 7.4 Fitting diagonal/off-diagonal branches to the lattice

- `_fit_to_lattice_diag`: for a diagonal branch, the MPO tensor is placed at the position on only one branch, and the `vd2 ⊗ I₂` identity tensor is placed at the position on the other branch;
- `_fit_to_lattice_offdiag`: for an off-diagonal branch, the `u` block is placed at `pos1` (= the minimum position) and the `v` block at `pos2`, with identity tensors filling the positions in between (`band_boundary` gives the position intervals of the ± branches for each time step).

In this way, the exponential gates of the off-diagonal coupling are "split" onto the lattice positions of the upper and lower branches, which is physically equivalent to the copy tensor of the effective-Hamiltonian path.

### 7.5 PT→ADT conversion

`copytensor` in `conversions.jl` constructs the copy tensor `m[i,i,i]=1`; the full PT→ADT conversion (`toadt`) is commented out/a placeholder in the current version (implementation hints for the method of the paper can be found in the comments at the head of the file).

## 8. System dynamics and boundary conditions

- `boundarycondition!` (`src/boundarycondition.jl`): applies the initial state (`initialstate!`) at the lattice boundaries. On the real-time contour, the Grassmann representation of the initial density matrix is placed at positions 1 and 2; in imaginary time, vacuum/thermal-state boundaries are placed at the two ends;
- `sysdynamics`/`sysdynamics!` (`src/models/`): multiplies the system time-evolution operators (first-order exponential gates `exp(-iH δt)`, or Lindblad propagators) into the ADT/PT, corresponding to the system propagators of Eqs. (9)-(10) of the paper;
- `ImpurityHamiltonian`/`ImpurityLindbladian`: defined in `ImpurityModelBase` as the system Hamiltonian/dissipative operators; `sysdynamics` generates the propagators from them and applies them in place.

## 9. The MPO Hamiltonian machinery

`MPOHamiltonian` (`src/mpohamiltonian/mpohamiltonian.jl`) is a column of `SchurMPOTensor`s, used for:

- building the MPO representation of the system Hamiltonian (`tompotensors`);
- time evolution `timeevompo(m, dt; alg=WII())` (site-by-site evolution);
- long-range coupling terms: `GenericDecayTerm` (generic exponential decay), `PowerlawDecayTerm` (power law, e.g., long-range XXZ).

`get_A/B/C/D` extracts the block structure of a `SchurMPOTensor`, and `_SiteW_impl` assembles it into a `SparseMPOTensor`.

## 10. Tensor-network algorithms

### 10.1 Orthogonalization

(`src/adt/orth.jl`) `Orthogonalize` is configured with `orth` (QR/SVD), `trunc`, `normalize`, `verbosity`:

- `leftorth!`: applies `tqr!`/`tsvd!` from left to right; after R/S normalization, the factor is absorbed into the next site;
- `rightorth!`: applies `tlq!`/`tsvd!` from right to left; after L normalization, the factor is absorbed into the previous site;
- `canonicalize!`: first left-orthogonalizes with QR (without truncation), then right-orthogonalizes with the SVD specified by `alg` and truncates. The comments explicitly warn: for ADT/ProcessTensor, **do not enable normalize** (renormalization would break the scaling semantics of the influence functional).

On the QR path, truncation has no effect (a `@warn` is issued). The truncation error is reported as the relative error `sqrt(ε²/(n²+ε²))`.

### 10.2 MPS multiplication `mult!` / `mult`

(`src/adt/mult/svdmult.jl`) the `mult!(x, y)` algorithm:

1. `tqr!` at the left end (`tie(n_fuse(...), ...)` merges the indices);
2. iterate from left to right: contract `tmp = r ⊗ x[i] ⊗ y[i]` → `n_fuse` → `tqr!` → `_renormalize!`;
3. after finishing at the right end, `_rightorth!(x, SVD(), trunc)` performs the SVD truncation from right to left;
4. `setscaling!(x, scaling(x)*scaling(y))`.

`mult(x, y) = mult!(copy(x), y)` is the non-mutating version. `DMRGMult1` (a `DMRGAlgorithm`) provides a variant with an `initguess` (default `:svd`), combined with the `D`/`tol` truncation of `SVDCompression`.

### 10.3 Truncation schemes

- `truncdim(D)`: truncates only by the maximum bond dimension;
- `trunccutoff(ε)`: truncates according to the accumulated truncation error;
- `truncdimcutoff(; D, ε)`: combines the two; `NoTruncation()` performs no truncation.

### 10.4 Transfer matrices

`ADTTransferMatrix` (`src/observables/adt/transfer.jl`) assembles a set of ADTs into a transfer operator:

- `left * m`: contracts `transfer_left` site by site from left to right (with `lmul!(scaling, …)` accounting for the scaling);
- `m * right`: from right to left;
- `TransferMatrix(states...)` / `TransferMatrix(j, states...)`: the whole-chain or single-site versions (used for single-site operators).

## 11. Observable evaluation

### 11.1 Environment caches

`environments(lattice, A, B...)` (`src/observables/`) constructs the environment caches:

- `ADTExpectationCache{A, Bs, lattice, hleft, hright}`: `hleft[i]` is the left environment of the first i sites (accumulated from left to right via `left * TransferMatrix(i, As...)`), and `hright[i]` is the right environment; `Zvalue = only(hleft[end])`;
- `PTExpectationCache`: the PT version; on the real-time contour it accepts the keyword `ρ₀` (initial density matrix) as the right boundary condition.

### 11.2 Expectation values

- `expectation(m, cache)`: applies the operator `m` (`apply!`) to `copy(cache.A)`, contracts the transfer matrices from right to left over the support of the operator, and finally contracts the two end environments with `contract_center(left, right)`;
- `expectationvalue(m, cache) = expectation(m, cache) / Zvalue(cache)`: the normalized expectation value $\langle m\rangle = \mathrm{tr}_{\text{bath}}[\bar\psi m \psi]/Z$;
- `Zvalue(cache)`: the partition function Z.

## 12. Error sources and hyperparameters

The five error sources of the paper and their counterparts in the code:

| Error source | Hyperparameters | Description |
|---|---|---|
| Time discretization | `δt` (`δτ`) | First-order Trotter / influence-functional discretization error |
| SVD truncation | `truncdim` / `trunccutoff` / `truncdimcutoff` | Bond dimension `χ` |
| Exponential expansion | `algexpan` | Number of terms and tolerance of the Prony expansion (default `OverDeterminedProny(n=15, tol=1e-4)`) |
| Translationally invariant refinement | `k` (default 5), `fast` | With `fast=true`: first build the differential IF of width `dt/2^k`, then square it k times by tree bisection to obtain the full-length influence functional |
| System propagator | `algevo` (`WII()`), `algmult` (`DefaultMultAlg`) | Accuracy of the MPO time evolution and of the multiplication compression |

The defaults are given in `src/defaults.jl`: `DefaultTruncation = truncdimcutoff(D=100, ϵ=1e-14)`, `DefaultITruncation = truncdimcutoff(D=200, ϵ=1e-10)`, `DefaultKTruncation = truncdimcutoff(D=1000, ϵ=1e-10)`, `DefaultIntegrationTruncation = DefaultMPOTruncation = truncdimcutoff(D=10000, ϵ=1e-12)`; `XTRGIF(; algexpan=OverDeterminedProny(n=15, tol=1e-4), algevo=WII(), algmult=DefaultMultAlg, k=5, fast=true)`.

## 13. Overview of the data flow

To compute a complete time-evolution dynamics:

```
corr = bath correlation function (ImpurityModelBase spectral function/correlation function)
  └─ exponential_expansion → ηᵢⱼ coefficients
lattice = ADTLattice/PTLattice (contour + Fock ordering + index mapping)
hyb = AdditiveHyb (diagonal) or NonDiagonalHyb/GeneralHybStyle (off-diagonal)
  ├─ diagonal: PartialIF (bond dim 2) row-by-row multiplication, or XTRGIF (MPO exponential expansion + WI/WII evolution)
  └─ off-diagonal: PT path, 4 branch MPOs (η⁺⁺…η⁻⁻) + fused_op/split_mpotensor + branch-wise multiplication
sysdynamics!: multiply in the system propagator (H or Lindblad)
boundarycondition!: initial state
observables: environments cache → expectationvalue/expectation/Zvalue
```

Scaling of the core loop: the number of sites `L` of the ADT/PT grows linearly with the number of time steps, and the bond dimension is controlled by the truncation scheme; the TTI scheme replicates the single-step differential influence functional (an MPO) onto all sites in a translationally invariant pattern, greatly reducing the cost of constructing the influence functional.
