# INTERVAL TYPES

"""
Abstract interval type
"""
abstract type AbstractInterval{T<:AbstractFloat} end

const PersistentDiagram{T} = AbstractVector{<:AbstractInterval{T}}

"""
    first(i::AbstractInterval)

Return a birth value of the interval `i`.
"""

first(i::AbstractInterval) = error("`first` is not implemented ")
"""
    last(i::AbstractInterval)

Return a death value of the interval `i`.
"""
last(i::AbstractInterval) = error("`last` is not implemented ")

"""
    last(i::AbstractInterval)

Return a death value of the interval `i`.
"""
dim(i::AbstractInterval) = error("`dim` is not implemented ")

# auxilary methods
show(io::IO, i::AbstractInterval) = print(io, "[$(first(i)),$(last(i)))")
isless(i1::AbstractInterval, i2::AbstractInterval) = first(i1) < first(i2) ? true : ( first(i1) == first(i2) ? last(i1) < last(i2) : false )
birth(i::AbstractInterval) = first(i) - last(i)
death(i::AbstractInterval) = first(i) + last(i)
pair(i::AbstractInterval) = first(i) => last(i)

"""
Simple implementation of the `AbstractInterval` type
"""
struct Interval{T<:AbstractFloat} <: AbstractInterval{T}
    dim::Int
    b::T
    d::T
end
Interval(dim::Int, p::Pair{T,T}) where {T<:AbstractFloat} = Interval(dim, first(p), last(p))
Interval(p::Pair) = Interval(0, p)
intervals(d::Int, ps::Pair...) = [Interval(d, p) for p in ps]

dim(i::Interval) = i.dim
first(i::Interval) = i.b
last(i::Interval) = i.d
diag(i::Interval) = let c = death(i)/2.0; Interval(dim(i), c, c) end

"""
Interval annotated with a generator
"""
struct AnnotatedInterval{T<:AbstractFloat} <: AbstractInterval{T}
    dim::Int
    b::T
    d::T
    generator::AbstractChain
end
AnnotatedInterval(dim::Int, b::T, d::T) where {T<:AbstractFloat} = AnnotatedInterval(dim, b, d, EmptyChain())
