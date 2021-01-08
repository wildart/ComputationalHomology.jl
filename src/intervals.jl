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
isempty(i::AbstractInterval) = birth(i) == death(i)
eltype(i::AbstractInterval{T}) where {T<:AbstractFloat} = T
in(e::T, i::AbstractInterval{T}) where {T<:AbstractFloat} = birth(i) <= e < death(i)
(==)(i1::AbstractInterval, i2::AbstractInterval) =
    birth(i1) == birth(i2) && death(i1) == death(i2)
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
show(io::IO, i::AnnotatedInterval) = print(io, "[$(i.b),$(i.d)) => $(i.g)")

"""
    diagram(d, pts)

Construct persistence diagram of dimension `d` from point pairs `pts`.
"""
diagram(pts::Pair...) = [Interval(p) for p in pts]

#
# I/O
#

function write(io::IO, pdd::Dict{Int, PD}) where {T<:Real, PD<:PersistenceDiagram{T}}
    for (d, pd) in pdd
        write(io, "$d\n")
        for i in pd
            write(io, "$(birth(i)) $(death(i))\n")
        end
    end
end

function read(io::IO, ::Type{I}) where {T<:Real, I<:AbstractInterval{T}}
    pd = Dict{Int, Vector{I}}()
    d = -1
    while !eof(io)
        l = readline(io)
        vals = split(l, ' ')
        if length(vals) == 1
            d = parse(Int, vals[1])
            pd[d] = I[]
        else
            vs = parse.(T, vals)
            push!(pd[d], Interval(vs[1],vs[2]))
        end
    end
    return pd
end
