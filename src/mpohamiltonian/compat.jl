# FiniteMPSAlgorithms 后端的兼容层
#
# 旧版 TEMPO 的 `TimeEvoMPOAlgorithm` 是本地抽象类型（WI/WII/ComplexStepper
# 的公共父类）；FiniteMPSAlgorithms 中这些 stepper 直接继承 `MPSAlgorithm`，
# 因此这里用类型别名保持 `XTRGIF` 等接口的类型约束不变。

const TimeEvoMPOAlgorithm = MPSAlgorithm

# timeevompo 返回 plain `MPO`（FMA 的 `tompotensors` 只收 `MPOHamiltonian`），
# 这里补上到稠密 4 指标站点张量列表的转换，保持与 `MPOHamiltonian` 版本一致
tompotensors(h::MPO) = h.data
