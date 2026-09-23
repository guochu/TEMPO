# FiniteMPSAlgorithms 后端的兼容层
#
# 旧版 TEMPO 的 `TimeEvoMPOAlgorithm` 是本地抽象类型（WI/WII/ComplexStepper
# 的公共父类）；FiniteMPSAlgorithms 中这些 stepper 直接继承 `MPSAlgorithm`，
# 因此这里用类型别名保持 `XTRGIF` 等接口的类型约束不变。

const TimeEvoMPOAlgorithm = MPSAlgorithm
