# MPO 哈密顿量与稀疏 MPO 张量（后端：FiniteMPSAlgorithms）
#
# AbstractSparseMPOTensor / SparseMPOTensor / SchurMPOTensor / MPOHamiltonian /
# tompotensors 以及 WI / WII / ComplexStepper / FirstOrderStepper / complex_stepper /
# timeevompo 全部来自 FiniteMPSAlgorithms（弱绑定直接使用并随 TEMPO re-export）。
# 本目录保留 TEMPO 特有的长程衰减项（ExponentialDecayTerm / GenericDecayTerm /
# PowerlawDecayTerm、expand_decayterm）以及少量兼容构造器（compat.jl）。

include("schurmpo/schurmpo.jl")
include("compat.jl")
