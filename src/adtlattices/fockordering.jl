"""
    FockOrdering

Abstract type for Fock space orderings, specifying the ordering convention of fermionic (Grassmann) operators on the lattice.

Concrete orderings are given by the subtypes: [`ImagFockOrdering`](@ref) (imaginary-time contour), [`RealFockOrdering`](@ref) (real-time contour), and [`MixedFockOrdering`](@ref) (mixed contour).
"""
abstract type FockOrdering end
"""
    ImagFockOrdering <: FockOrdering

Abstract Fock ordering on the imaginary-time (Matsubara) contour.
"""
abstract type ImagFockOrdering <: FockOrdering end
"""
    RealFockOrdering <: FockOrdering

Abstract Fock ordering on the real-time (Keldysh) contour.
"""
abstract type RealFockOrdering <: FockOrdering end
"""
    MixedFockOrdering <: FockOrdering

Abstract Fock ordering on the mixed (Kadanoff-Baym, real-time + imaginary-time) contour.
"""
abstract type MixedFockOrdering <: FockOrdering end


"""
    LayoutStyle

Abstract type for lattice layout styles, specifying how states are arranged in the tensor network.
"""
abstract type LayoutStyle end
"""
    TimeLocalLayout <: LayoutStyle

Time-local layout: the state at each time step remains locally arranged, suitable for acting with time-local operators.
"""
struct TimeLocalLayout <: LayoutStyle end


"""
    TimeOrderingStyle

Abstract type for time-ordering styles, specifying the time-ordering direction on the real/imaginary-time branches.

- `TimeAscending`: ascending time ordering;
- `TimeDscending`: descending time ordering.
"""
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

First-order Fock ordering on the imaginary-time contour: each imaginary-time step contains two orders (M2, M1),
with descending time ordering (`TimeDscending`) on the imaginary-time branch. Alias `MM`.
"""
struct M2M1 <: ImagFockOrdering end
LayoutStyle(::Type{M2M1}) = TimeLocalLayout()
ImaginaryTimeOrderingStyle(::Type{<:ImagFockOrdering}) = TimeDscending()
TimeOrderingStyle(::Type{O}) where {O<:ImagFockOrdering} = ImaginaryTimeOrderingStyle(O)
const MM = M2M1

"""
	struct M2m2M1m1 <: RealFockOrdering 

First-order Fock ordering on the real-time contour: each real-time step contains two orders (M2, M1) and their conjugates (m2, m1),
with descending time ordering (`TimeDscending`) on the real-time branch. Alias `MmMm`.
"""
struct M2m2M1m1 <: RealFockOrdering end
LayoutStyle(::Type{M2m2M1m1}) = TimeLocalLayout()
RealTimeOrderingStyle(::Type{<:RealFockOrdering}) = TimeDscending()
TimeOrderingStyle(::Type{O}) where {O<:RealFockOrdering} = RealTimeOrderingStyle(O)
const MmMm = M2m2M1m1

"""
	struct M2M1_m1M1m2M2 <: MixedFockOrdering

Fock ordering on the mixed contour: ascending time ordering (`TimeAscending`) on the real-time branch and
descending time ordering (`TimeDscending`) on the imaginary-time branch.
"""
struct M2M1_m1M1m2M2 <: MixedFockOrdering end
LayoutStyle(::Type{M2M1_m1M1m2M2}) = TimeLocalLayout()
RealTimeOrderingStyle(::Type{M2M1_m1M1m2M2}) = TimeAscending()
ImaginaryTimeOrderingStyle(::Type{<:MixedFockOrdering}) = TimeDscending()
