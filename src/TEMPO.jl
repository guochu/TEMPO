module TEMPO

# auxiliary — tensor primitives（后端：FiniteMPSAlgorithms）
export TruncationScheme, NoTruncation, TruncateDimCutoff, truncdimcutoff, truncdim, truncrelerr, renyi_entropy
export TruncateDim, TruncateRelError
export SVDCompression
export OrthogonalFactorizationAlgorithm, leftorth!, rightorth!, leftorth, rightorth, tsvd!, tsvd, permute, isometry
# ContourIndex
export ContourIndex, branch, scalartype
# MPOHamiltonian
export MPOHamiltonian, tompotensors, timeevompo, WI, WII, ComplexStepper, FirstOrderStepper, complex_stepper
export SchurMPOTensor, SparseMPOTensor, ExponentialDecayTerm, GenericDecayTerm, PowerlawDecayTerm
export expand_decayterm

# ADT
export space_l, space_r, bond_dimension, bond_dimensions, scaling, phydim, phydims
export ADT, randomadt, isleftcanonical, isrightcanonical, iscanonical
export distance, distance2, Orthogonalize, leftorth!, rightorth!, canonicalize!
export mult, mult!, DMRG1
# PT
export ProcessTensor, randompt
# ADT terms
export ADTTerm, apply!
export AbstractFockTerm, FockTermS, FockTerm, ProdFockTerm
# ADT Lattices
export FockOrdering, ImagFockOrdering, RealFockOrdering, MixedFockOrdering, TimeOrderingStyle
export M2M1, M2m2M1m1, M2M1_m1M1m2M2, TimeLocalLayout
export index, OrderingStyle, LayoutStyle, ImaginaryTimeOrderingStyle, RealTimeOrderingStyle
export AbstractADTLattice, branches, ADTLattice, vacuumstate, indexmappings
export ImagADTLattice, ImagADTLattice1Order, RealADTLattice, RealADTLattice1Order
export MixedADTLattice, MixedADTLattice1Order
# PT Lattices
export AbstractPTLattice, PTLattice, integrate, ContourOperator, correlationfunction
export ImagPTLattice, ImagPTLattice1Order, RealPTLattice, RealPTLattice1Order
export MixedPTLattice, MixedPTLattice1Order
export rdm, quantummap, meanforcestate, mfs
# influence functional
export HybridizationStyle, AdditiveHyb, NonAdditiveHyb, NonDiagonalHyb, pairop
export PartialIF, XTRGIF, TDVPIF
export influenceoperators, influenceoperatorsteppers, influenceoperatorstepper
export hybriddynamics, hybriddynamics!, hybriddynamics_naive, hybriddynamics_naive!
export partialif_naive, partialif
# boundary condition
export boundarycondition, boundarycondition!, initialstate!
# models
export sysdynamics, sysdynamics!, ImpurityHamiltonian, QuenchedImpurityHamiltonian, TdImpurityOp, TdImpurityHamiltonian, ImpurityLindbladian
# observables
export environments, expectationvalue, Zvalue, expectation, TransferMatrix
export l_LL, r_RR




using Reexport
@reexport using ExpExp

using ImpurityModelBase, QuAPI
using KrylovKit: exponentiate, Arnoldi
import QuAPI: branch, index
using LinearAlgebra
using Base: @boundscheck
using TensorOperations,TupleTools
const TO = TensorOperations

# ---------------------------------------------------------------------------
# MPS/MPO 运算后端：FiniteMPSAlgorithms
#
# 张量运算原语（截断方案、tsvd!/leftorth!、tie/permute/isometry…）、
# CanonicalMPS/CanonicalMPO 链类型、ALS 迭代引擎（MultCache/HadamardCache）、
# SchurMPOTensor/SparseMPOTensor/MPOHamiltonian 以及 WI/WII/ComplexStepper
# 演化全部来自 FiniteMPSAlgorithms；TEMPO 的 ADT/ProcessTensor 等公共接口
# 是其实现之上的轻量级 wrapper（见 fmabackend.jl 与 adt/ pt/ mpohamiltonian/）。
#
# `mult`、`canonicalize!`、`scaling` 等名字通过 `import` 扩展（单一函数，
# TEMPO 方法与 FiniteMPSAlgorithms 方法共存于同一泛型函数）；
# `DMRG1`、`Defaults` 等保留 TEMPO 本地定义（遮蔽 FiniteMPSAlgorithms 的弱绑定）。
# ---------------------------------------------------------------------------
using FiniteMPSAlgorithms
# TEMPO 扩展（定义了新方法）的函数：必须 import，与 FMA 共用同一泛型函数
import FiniteMPSAlgorithms: mult, mult!, canonicalize!, canonicalize, leftorth!, rightorth!,
	swap!, permute!, permute, distance, distance2, scaling, setscaling!,
	svectors_uninitialized, unset_svectors!, changebond!, phydim, SVDCompression
# 仅调用（无 TEMPO 方法）的内部原语与类型：纯引入
using FiniteMPSAlgorithms: _renormalize!, _reduce_hadamard_site, _reduce_site, _env_updateright,
	_updateright, _contract_last,
	MPSAlgorithm, SchurMPOTensor, SparseMPOTensor, MPOHamiltonian,
	QR, QRpos, LQ, LQpos, SVD, SDD, Polar,
	OrthogonalFactorizationAlgorithm

# mps algorithms（TEMPO 公共算法配置；DMRG1 兼容层等）
include("algorithms.jl")

# default constants
include("defaults.jl")


# mpo hamiltonian
include("mpohamiltonian/mpohamiltonian.jl")


include("contourindices.jl")

# adt
include("adt/adt.jl")

# pt
include("pt/pt.jl")

# FiniteMPSAlgorithms 适配层（零拷贝视图、DMRG1 翻译等；依赖 ADT/ProcessTensor）
include("fmabackend.jl")

# conversion between pt and adt
include("conversions.jl")


# adtterms
include("adtterms.jl")
include("fockterms.jl")


# adtlattices
include("adtlattices/adtlattices.jl")

# ptlattices
include("ptlattices/ptlattices.jl")

include("contouroperators.jl")

# correlation function
include("correlationfunction.jl")

# Feynman-Vernon influence functional
include("influencefunctional/influencefunctional.jl")

# Time-dependent Feynman-Vernon influence functional
include("tdinfluencefunctional/tdinfluencefunctional.jl")


# boundary condition
include("boundarycondition.jl")

# models
include("models/models.jl")

# observables
include("observables/observables.jl")

end
