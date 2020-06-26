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

reduce(::Type{R}, ∂::Vector{<:AbstractSet}) where {R<:AbstractPersistenceReduction} = reduce!(R, deepcopy(∂))

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
function diagram(flt::Filtration{C,FI}, R::Vector, length0=false,
                 absolute=true) where {C <: AbstractComplex, FI <: AbstractFloat}
    # resulting intervals
    intrs = Dict{Int,Vector{Interval{FI}}}()

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
                    push!(intrs[sdim], Interval(bv, dv))
                end
            end
        end
    end
    for i in births
        sdim, si, fv = ord[i]
        !haskey(intrs, sdim) && setindex!(intrs, Interval[], sdim)
        push!(intrs[sdim], Interval(fv, Inf))
    end

    return intrs
end

diagram(::Type{R}, flt::Filtration; length0=false, absolute=true) where {R <: AbstractPersistenceReduction} =
    diagram(flt, reduce!(R, boundary(flt)), length0, absolute)
diagram(flt::Filtration; kwargs...) = diagram(TwistReduction, flt; kwargs...)

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
diagram(h::PersistentHomology{R}) where {R <: AbstractPersistenceReduction} = diagram(R, h.filtration)

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


#
# Persistent Cohomology
#

function coboundarymap(coch::Chain{IX,R}, ∂::Chain{IX,R}) where {IX<:Integer, R}
    sum(coch[be]*bc for (be,bc) in ∂)
end

function cocycles(ccycs, ∂::Chain{IX,R}, w::Float64) where {IX<:Integer, R}
    cᵢ = zeros(R, length(ccycs))
    for (i,intr) in enumerate(ccycs)
        cᵢ[i] = if w ∉ intr
            zero(R)
        else
            coboundarymap(intr.g, ∂)
        end
    end
    cᵢ
end

function update!(cyc::Dict{Int,Vector{AnnotatedInterval{F,Chain{IX,R}}}},
                 cplx, splxs, w::F) where {F<:AbstractFloat, IX<:Integer, R}
    d = first(first(splxs))
    d >= length(cyc) && return
    CXₗ₋₁ = cyc[d-1]
    for splx in splxs
        σ = cplx[last(splx), d]
        dαᵢ = cocycles(CXₗ₋₁, boundary(R, σ), w)
        idxs = findall(!iszero, dαᵢ)
        if !isempty(idxs)
            j = idxs[end]
            cⱼ = dαᵢ[j]
            αⱼ = CXₗ₋₁[j].g
            for i in idxs[1:end-1]
                cᵢ = dαᵢ[i]
                CXₗ₋₁[i].g -= (cᵢ/cⱼ)*αⱼ
            end
            CXₗ₋₁[j].d = w
        else
            intr = AnnotatedInterval(w, F(Inf), Chain(d, [splx[2]], [one(R)]))
            push!(cyc[d], intr)
        end
    end
end

function persistentcohomology(::Type{R}, flt::Filtration;
                              length0=false, maxoutdim=dim(complex(flt))) where {R}
    cplx = complex(flt)
    d = dim(cplx)
    itr = Iterators.flatten(
            ((v => filter(e->first(e) == i, splxs))
                for i in (:)(extrema(map(first, splxs))...)) for (v, splxs) in flt )

    cyc = Dict(i=>AnnotatedInterval{R,Chain{UInt64,R}}[] for i in -1:min(d-1, maxoutdim))
    for (w, splxs) in itr
        update!(cyc, cplx, splxs, w)
    end
    delete!(cyc, -1)

    # remove all intervals of 0-length
    cyc = !length0 ? Dict( d=>filter!(!isempty, c) for (d,c) in cyc) : cyc

    return cyc
end
persistentcohomology(flt::Filtration; kwargs...) = persistentcohomology(Float64, flt; kwargs...)
