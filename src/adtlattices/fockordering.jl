abstract type FockOrdering end
abstract type ImagFockOrdering <: FockOrdering end
abstract type RealFockOrdering <: FockOrdering end
abstract type MixedFockOrdering <: FockOrdering end


abstract type LayoutStyle end
struct TimeLocalLayout <: LayoutStyle end


abstract type TimeOrderingStyle end
struct TimeAscending <: TimeOrderingStyle end
struct TimeDscending <: TimeOrderingStyle end

LayoutStyle(x::FockOrdering) = LayoutStyle(typeof(x))
# similargrassmannordering(::Type{T}) where {T<:FockOrdering} = error("similargrassmannordering not implemented for FockOrdering type $T")
# similargrassmannordering(x::FockOrdering) = similargrassmannordering(typeof(x))

ImaginaryTimeOrderingStyle(x::FockOrdering) = ImaginaryTimeOrderingStyle(typeof(x))
RealTimeOrderingStyle(x::FockOrdering) = RealTimeOrderingStyle(typeof(x))
TimeOrderingStyle(x::FockOrdering) = TimeOrderingStyle(typeof(x))

"""
	struct M2M1 <: ImagFockOrdering
"""
struct M2M1 <: ImagFockOrdering end
LayoutStyle(::Type{M2M1}) = TimeLocalLayout()
ImaginaryTimeOrderingStyle(::Type{<:ImagFockOrdering}) = TimeDscending()
TimeOrderingStyle(::Type{O}) where {O<:ImagFockOrdering} = ImaginaryTimeOrderingStyle(O)
const MM = M2M1

"""
	struct M2m2M1m1 <: RealFockOrdering 
"""
struct M2m2M1m1 <: RealFockOrdering end
LayoutStyle(::Type{M2m2M1m1}) = TimeLocalLayout()
RealTimeOrderingStyle(::Type{<:RealFockOrdering}) = TimeDscending()
TimeOrderingStyle(::Type{O}) where {O<:RealFockOrdering} = RealTimeOrderingStyle(O)
const MmMm = M2m2M1m1

"""
	struct M2M1_m1M1m2M2 <: MixedFockOrdering
"""
struct M2M1_m1M1m2M2 <: MixedFockOrdering end
LayoutStyle(::Type{M2M1_m1M1m2M2}) = TimeLocalLayout()
RealTimeOrderingStyle(::Type{M2M1_m1M1m2M2}) = TimeAscending()
ImaginaryTimeOrderingStyle(::Type{<:MixedFockOrdering}) = TimeDscending()
