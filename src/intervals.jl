# INTERVAL TYPES

"""
Abstract interval type
"""
abstract type AbstractInterval{T<:AbstractFloat} end

const PersistenceDiagram{T} = AbstractVector{<:AbstractInterval{T}}

"""
    birth(i::AbstractInterval)

Return a birth value of the interval `i`.
"""
birth(i::AbstractInterval) = error("`birth` is not implemented ")

"""
    death(i::AbstractInterval)

Return a death value of the interval `i`.
"""
death(i::AbstractInterval) = error("`death` is not implemented ")

# auxilary methods
show(io::IO, i::AbstractInterval) = print(io, "[$(birth(i)),$(death(i)))")
birthx(i::AbstractInterval) = birth(i) - death(i)
deathx(i::AbstractInterval) = birth(i) + death(i)
pair(i::AbstractInterval) = birth(i) => death(i)
in(e::T, i::AbstractInterval{T}) where {T<:AbstractFloat} = birth(i) <= e < death(i)
isless(i1::AbstractInterval, i2::AbstractInterval) =
    birth(i1) < birth(i2) ? true : ( birth(i1) == birth(i2) ? death(i1) < death(i2) : false )

"""
    Interval{T<:AbstractFloat} <: AbstractInterval{T}

Simple implementation of the `AbstractInterval` type.
"""
struct Interval{T<:AbstractFloat} <: AbstractInterval{T}
    b::T
    d::T
end
Interval(p::Pair{T,T}) where {T<:AbstractFloat} = Interval(p.first, p.second)
Interval() = Interval(0.0=>Inf)

birth(i::Interval) = i.b
death(i::Interval) = i.d
diag(i::Interval) = let c = deathx(i)/2.0; Interval(c, c) end

"""
    AnnotatedInterval{T<:AbstractFloat, C<:AbstractChain} <: AbstractInterval{T}

Interval type annotated with a generator chain.
"""
mutable struct AnnotatedInterval{T<:AbstractFloat, C<:AbstractChain} <: AbstractInterval{T}
    b::T
    d::T
    g::C
end
AnnotatedInterval(b::T, d::T) where {T<:AbstractFloat} = AnnotatedInterval(b, d, EmptyChain())

birth(i::AnnotatedInterval) = i.b
death(i::AnnotatedInterval) = i.d
generator(i::AnnotatedInterval) = i.g

"""
    diagram(d, pts)

Construct persistence diagram of dimension `d` from point pairs `pts`.
"""
diagram(pts::Pair...) = [Interval(p) for p in pts]
