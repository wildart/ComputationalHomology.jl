import LinearAlgebra: diag

# INTERVAL TYPES

"""
Abstract interval type
"""
abstract type AbstractInterval{T<:AbstractFloat} end

const PersistentDiagram{T} = AbstractVector{AbstractInterval{T}}

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

# REDUCTION ALGORITHMS

abstract type AbstractPersistenceReduction end
struct StandardReduction <: AbstractPersistenceReduction end
struct TwistReduction <: AbstractPersistenceReduction end

# Boundary matrix reduction algorithms
lastindex(col::AbstractSet) = length(col) == 0 ? -1 : maximum(col)

"""Standart reduction"""
function reduce!(::Type{StandardReduction}, ∂::Vector{<:AbstractSet})
    lowest_one_lookup = fill(-1, length(∂))
    for col in eachindex(∂)
        lowest_one = lastindex(∂[col])
        while lowest_one != -1 && lowest_one_lookup[lowest_one] != -1
            lowest = lowest_one_lookup[lowest_one]
            symdiff!(∂[col], ∂[lowest])
            lowest_one = lastindex(∂[col])
        end
        if lowest_one != -1
            lowest_one_lookup[lowest_one] = col
        end
    end
    return ∂
end

"""Twist reduction"""
function reduce!(::Type{TwistReduction}, ∂::Vector{<:AbstractSet})
    lowest_one_lookup = fill(-1, length(∂))

    for dim in maximum(map(length, ∂)):-1:1
        for col in eachindex(∂)
            if length(∂[col]) == dim
                lowest_one = lastindex(∂[col])
                while lowest_one != -1 && lowest_one_lookup[lowest_one] != -1
                    lowest = lowest_one_lookup[lowest_one]
                    symdiff!(∂[col], ∂[lowest])
                    lowest_one = lastindex(∂[col])
                end
                if lowest_one != -1
                    lowest_one_lookup[lowest_one] = col
                    empty!(∂[lowest_one])
                end
            end
        end
    end
    return ∂
end

Base.reduce(::Type{R}, ∂::Vector{<:AbstractSet}) where {R<:AbstractPersistenceReduction} = reduce!(R, deepcopy(∂))

function generate_pairs(∂::Vector{<:AbstractSet}; reduced = false)
    ridx = reduced ? 1 : 0
    births = eltype(∂)()
    ps = Pair[]
    for i in eachindex(∂)
        if length(∂[i]) > 0
            b = last(∂[i])
            d = i
            delete!(births, b)
            delete!(births, d)
            (d > b) && push!(ps, (b-ridx) => (d-ridx) )
        else
            push!(births, i)
        end
    end
    for i in births # no lowest, create semi-infinite interval
        push!(ps, (i-ridx)=>Inf)
    end
    return ps
end

"Compute raw persistence pairs (boundary matrix is reduced in a process)"
function pairs(::Type{R}, ∂::Vector{<:AbstractSet}; reduced = false) where {R <: AbstractPersistenceReduction}
    reduce!(R, ∂) # reduce  boundary matrix
    return generate_pairs(∂, reduced=reduced), ∂  # generate pairs
end
pairs(::Type{R}, flt::Filtration; reduced = false) where {R <: AbstractPersistenceReduction} =
    pairs(R, boundary(flt, reduced = reduced), reduced = reduced)

"""Return persistent diagram (birth-death pairs) per dimension."""
function intervals(flt::Filtration, R::Vector, length0=false, absolute=true)
    # resulting intervals
    intrs = Dict{Int,Vector{Interval}}()

    # compute intervals
    births = Set{Int}()
    ord = order(flt)
    for (i, (sdim, si, fv)) in enumerate(ord)
        col = R[i]
        if length(col) == 0
            push!(births, i)
        else
            b = lastindex(col)
            d = i
            delete!(births, b)
            delete!(births, d)
            if d > b
                sdim = absolute ? ord[b][1] : ord[d][1]
                bv = ord[b][3]
                dv = ord[d][3]
                if dv > bv || length0
                    !haskey(intrs, sdim) && setindex!(intrs, Interval[], sdim)
                    push!(intrs[sdim], Interval(sdim, bv, dv))
                end
            end
        end
    end
    for i in births
        sdim, si, fv = ord[i]
        !haskey(intrs, sdim) && setindex!(intrs, Interval[], sdim)
        push!(intrs[sdim], Interval(sdim, fv, Inf))
    end

    return intrs
end

intervals(::Type{R}, flt::Filtration; length0=false, absolute=true) where {R <: AbstractPersistenceReduction} =
    intervals(flt, reduce!(R, boundary(flt)), length0, absolute)
intervals(flt::Filtration; kwargs...) = intervals(TwistReduction, flt; kwargs...)

"Calculate persistent Betti numbers for a filtration complex of dimension `dim`"
function betti(flt::Filtration, R::Vector, p::Int)
    maxdim =  dim(complex(flt))
    @assert maxdim >= p "Cannot calculate $p-dimensional Betti number for $maxdim-complex"

    # indexes of p-simplices
    p2i = findall(d->d[1]==p, order(flt))

    # the number of zero columns that correspond to p-simplices
    z = sum(map(length, R[p2i]) .== 0)

    # the number of lowest ones in rows that correspond to p-simplices
    l = length( findall(l->l ∈ p2i, map(lastindex, R)) )

    β = z - l
    return β < 0 ? 0 : β
end
betti(::Type{R}, flt::Filtration, p::Int) where {R <: AbstractPersistenceReduction} =
    betti(flt, reduce!(R, boundary(flt)), p)
betti(flt::Filtration, p::Int) = betti(TwistReduction, flt, p)

"""
Persistent homology group iterator for a filtration
"""
mutable struct PersistentHomology{R <: AbstractPersistenceReduction} <: AbstractHomology
    filtration::Filtration
    ∂::Vector{<:AbstractSet}
end
function persistenthomology(::Type{R}, flt::Filtration;
                            reduced::Bool=false) where {R <: AbstractPersistenceReduction}
    RM = reduce!(R, boundary(flt, reduced = reduced))
    return PersistentHomology{R}(flt, RM)
end
persistenthomology(flt::Filtration) = persistenthomology(TwistReduction, flt)

show(io::IO, h::PersistentHomology{R}) where {R <: AbstractPersistenceReduction} =
    print(io, "PersistentHomology[$(h.filtration) with $R]")

#
# Interface methods
#

group(h::PersistentHomology, p::Int) = betti(h.filtration, h.∂, p)
intervals(h::PersistentHomology{R}) where {R <: AbstractPersistenceReduction} = intervals(R, h.filtration)

#
# Iterator methods
#

length(h::PersistentHomology) = dim(h.filtration.complex)+1

"""Return homology group type: dimension & Betti numbers."""
eltype(h::PersistentHomology) = Tuple{Int, Int}

function iterate(h::PersistentHomology, p=0)
    p > dim(h.filtration.complex) && return nothing
    βₚ = group(h, p)
    return (p, βₚ), p+1
end
