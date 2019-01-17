import LinearAlgebra: diag

abstract type AbstractPersistenceReduction end
mutable struct StandardReduction <: AbstractPersistenceReduction end
mutable struct TwistReduction <: AbstractPersistenceReduction end

struct Interval
    dim::Int
    generator::AbstractChain
    b::Number
    d::Number
end
Interval(dim::Int, b::Number, d::Number) = Interval(dim, EmptyChain(), b, d)
Interval(dim::Int, p::Pair) = Interval(dim, p.first, p.second)
Interval(p::Pair) = Interval(0, p)
Base.show(io::IO, intr::Interval) = print(io, "[$(intr.b),$(intr.d))")
Base.isless(i1::Interval, i2::Interval) = i1.b < i2.b ? true : ( i1.b == i2.b ? i1.d < i2.d : false )
birth(i::Interval) = i.b - i.d
death(i::Interval) = i.b + i.d
diag(i::Interval) = let c = death(i)/2.0; Interval(i.dim, i.generator, c, c) end

pair(i::Interval) = i.b => i.d
intervals(d::Int, ps::Pair...) = [Interval(d, p) for p in ps]

# Boundary matrix reduction algorithms
lastindex(col::BitSet) = length(col) == 0 ? -1 : last(col)

"""Standart reduction"""
function reduce!(::Type{StandardReduction}, ∂::Vector{BitSet})
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
function reduce!(::Type{TwistReduction}, ∂::Vector{BitSet})
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
                end
            end
        end
    end
    return ∂
end

Base.reduce(::Type{R}, ∂::Vector{BitSet}) where {R<:AbstractPersistenceReduction} = reduce!(R, deepcopy(∂))

function generate_pairs(∂::Vector{BitSet}; reduced = false)
    ridx = reduced ? 1 : 0
    births = BitSet()
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
function pairs(::Type{R}, ∂::Vector{BitSet}; reduced = false) where {R <: AbstractPersistenceReduction}
    reduce!(R, ∂) # reduce  boundary matrix
    return generate_pairs(∂, reduced=reduced), ∂  # generate pairs
end
pairs(::Type{R}, flt::Filtration; reduced = false) where {R <: AbstractPersistenceReduction} =
    pairs(R, boundary_matrix(flt, reduced = reduced), reduced = reduced)

"""Return persistent diagram (birth-death pairs) per dimension."""
function intervals(flt::Filtration, R::Vector, length0=false, absolute=true)
    # resulting intervals
    intrs = Dict{Int,Vector{Interval}}()

    # compute intervals
    births = BitSet()
    ord = order(flt)
    for (i, (sdim, si, fv)) in enumerate(ord)
        col = R[i]
        if length(col) == 0
            push!(births, i)
        else
            b = last(col)
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

intervals(flt::Filtration; reduction=TwistReduction, length0=false, absolute=true) =
    intervals(flt, reduce!(reduction, boundary_matrix(flt)), length0, absolute)

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
betti(flt::Filtration, p::Int; reduction=TwistReduction) =
    betti(flt, reduce!(reduction, boundary_matrix(flt)), p)

"""
Persistent homology group iterator for a filtration
"""
mutable struct PersistentHomology <: AbstractHomology
    filtration::Filtration
    reduction::DataType
    R::Vector{BitSet}
end
function persistenthomology(::Type{R}, flt::Filtration;
                            reduced::Bool=false) where R<:AbstractPersistenceReduction
    RM = reduce!(R, boundary_matrix(flt, reduced = reduced))
    return PersistentHomology(flt, R, RM)
end
persistenthomology(flt::Filtration) = persistenthomology(TwistReduction, flt)

Base.show(io::IO, h::PersistentHomology) = print(io, "PersistentHomology[$(h.filtration) with $(h.reduction)]")

"""Return homology group type: dimension & Betti numbers."""
Base.eltype(::Type{PersistentHomology}) = Tuple{Int, Int}

#
# Interface methods
#

group(h::PersistentHomology, p::Int) = betti(h.filtration, h.R, p)
intervals(h::PersistentHomology) = intervals(h.filtration, reduction=h.reduction)

#
# Iterator methods
#

Base.length(h::PersistentHomology) = dim(h.filtration.complex)+1
Base.eltype(h::PersistentHomology) = Tuple{Int, Int}
function Base.iterate(h::PersistentHomology, p=0)
    p > dim(h.filtration.complex) && return nothing
    βₚ = group(h, p)
    return (p, βₚ), p+1
end
