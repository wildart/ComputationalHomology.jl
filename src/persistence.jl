abstract type AbstractPersistenceReduction end
mutable struct StandardReduction <: AbstractPersistenceReduction end
mutable struct TwistReduction <: AbstractPersistenceReduction end

lastindex(col::IntSet) = length(col) == 0 ? -1 : last(col)

"""Standart reduction"""
function Base.reduce(::Type{StandardReduction}, ∂::Vector{IntSet})
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
function Base.reduce(::Type{TwistReduction}, ∂::Vector{IntSet})
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

function generate_pairs(∂::Vector)
    pairs = Pair{Int,Int}[]
    for col in eachindex(∂)
        if length(∂[col]) > 0
            birth = last(∂[col])
            death = col
            push!(pairs, birth=>death)
        end
    end
    return pairs
end

"Compute raw persistence pairs (boundary matrix is reduced in a process)"
function pairs(::Type{R}, ∂::Vector) where {R <: AbstractPersistenceReduction}
    reduce(R, ∂) # reduce  boundary matrix
    return generate_pairs(∂), ∂  # generate pairs
end

"""Return birth-death pairs per dimension"""
function intervals(flt::Filtration, ps::Vector{Pair{Int,Int}}; length0=false)
    ITT = keytype(flt.index)

    # Construct total order index of simplexes in filtration
    total = Dict{Int,Tuple{Int,Int,ITT}}() # idx => (dim, simplex_id, fvalue)
    idx = 1
    for fltval in sort!(collect(keys(flt.index)))
        for (d, ci) in sort(flt.index[fltval], lt=(x,y)->(x[1] < y[1]))
            total[idx] = (d, ci, fltval)
            idx+=1
        end
    end

    # construct intervals
    intrs = Dict{Int,Vector{Pair{ITT,ITT}}}()
    for (b,d) in ps
        (total[b][3] == total[d][3] && !length0) && continue # do not include 0-length intervals
        intr = total[b][3] => total[d][3]

        idim = total[d][1]-1
        !haskey(intrs, idim) && setindex!(intrs, Pair{ITT,ITT}[], idim)
        push!(intrs[idim], intr)
    end

    return intrs, total
end

function intervals(flt::Filtration; reduction=TwistReduction, length0=false, reduced = false)
    ps, _ = pairs(reduction, boundary_matrix(flt, reduced=reduced))
    return intervals(flt, ps, length0=length0)
end

"Calculate persistent Betti numbers for a filtration complex of dimension `dim`"
function betti(∂::Vector, R::Vector, p::Int)
    simdim(s) = length(s) == 0 ? 0 : length(s)-1
    sdims = map(simdim, ∂)
    @assert maximum(sdims) >= p "Cannot calculate $p-dimensional Betti number for $(maximum(sdims))-complex"

    # the number of zero columns that correspond to p-simplices
    # z = mapreduce(i->length(R[i]) == 0 ? 1 : 0, +, find(d->d==p, sdims))

    RR = spzeros(Int, length(R), length(R))
    for (j, I) in enumerate(R)
        lasti = 0
        for i in I
            lasti = i+1
            RR[lasti,j] = 1
        end
        if lasti > 0
            RR[lasti,j] = 2
        end
    end
    p2i = find(d->d==p, sdims)
    # the number of zero columns that correspond to p-simplices
    z = find(i->i==0, sum(RR[:,p2i],1)) |> length
    # the number of lowest ones in rows that correspond to p-simplices
    l = find(i->i==2, RR[p2i,:]) |> length
    β = z - l
    return β < 0 ? 0 : β
end

"""
Persistent homology group iterator for a filtration
"""
mutable struct PersistentHomology{G} <: AbstractHomology{G}
    filtration::Filtration
    reduction::DataType
    ∂::Vector{IntSet}
    R::Vector{IntSet}
end
persistenthomology{R <: AbstractPersistenceReduction, G}(flt::Filtration, ::Type{R}, ::Type{G}) = PersistentHomology{G}(flt, R, Vector{IntSet}(), Vector{IntSet}())
persistenthomology{R <: AbstractPersistenceReduction}(flt::Filtration, ::Type{R}) = persistenthomology(flt, R, Int)

Base.show(io::IO, h::PersistentHomology) = print(io, "PersistentHomology[$(h.filtration) with $(h.reduction)]")

"""Return homology group type: dimension & Betti numbers."""
Base.eltype{G}(::Type{PersistentHomology{G}}) = Tuple{Int, Int}

#
# Interface methods
#

function group{G}(h::PersistentHomology{G}, p::Int)
    if length(h.∂) == 0
        h.∂ = boundary_matrix(h.filtration)
    end
    if length(h.R) == 0
        h.R = reduce(h.reduction, deepcopy(h.∂))
    end
    return betti(h.∂, h.R, p)
end

#
# Iterator methods
#

Base.length{G}(h::PersistentHomology{G}) = dim(h.filtration.complex)+1
Base.start{G}(h::PersistentHomology{G}) = 0
Base.done{G}(h::PersistentHomology{G}, state) = dim(h.filtration.complex) < state[1]
function Base.next{G}(h::PersistentHomology{G}, state)
    p = state[1]
    βₚ = group(h, p)
    return (p, βₚ), p+1
end
