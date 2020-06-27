"""
    iterate(cplx::AbstractComplex) -> Union{Nothing, AbstractCell}

Advance the iterator to obtain the next cell in the complex `cplx`
"""
function iterate(cplx::AbstractComplex, (tid, d, did)=(0, -1, 1))
    tid >= length(cplx) && return nothing
    if d < 0
        d = 0
    end
    if did > size(cplx, d)
        d += 1
        did = 1
    end
    return cells(cplx, d)[did], (tid+1, d, did+1)
end

"""
    iterate(flt::Filtration{C,FI,IX}) -> Union{Nothing, Tuple{FI,Vector{Tuple{Int,IX}}}}

Advance the filtration iterator `flt` to obtain the next filtration value of type `FI` with
a corresponding collection of cell parameters: dimension & identifier of type `IX`.
"""
function iterate(flt::Filtration{C,FI,IX}, state=nothing) where {C<:AbstractComplex, FI<:AbstractFloat, IX<:Integer}
    ord = order(flt)
    if state === nothing # calculate initial state
        idx = 1
        fval = ord[idx][3]
        incr = (maximum(flt)-minimum(flt)) / flt.divisions
    else
        idx, fval, incr = state
    end
    idx > length(ord) && return nothing # done
    splxs = Tuple{Int,IX}[] #simplex dim & index
    while idx <= length(ord) && (fval+incr) >= ord[idx][3]
        push!(splxs, ord[idx][1:2])
        idx += 1
    end
    nextfval = fval+incr
    if idx <= length(ord) && isinf(flt.divisions)
        nextfval = ord[idx][3]
    end
    return (fval, splxs), (idx, nextfval, incr)
end

"""
Persistent homology group iterator for a filtration
"""
mutable struct PersistentHomology{R <: AbstractPersistenceReduction, C, FI, IX} <: AbstractHomology
    filtration::Filtration{C,FI,IX}
    ∂::Vector{<:AbstractSet}
    length0::Bool
end
function persistenthomology(::Type{R}, flt::Filtration{C,FI,IX}; length0=false,
                            reduced::Bool=false) where {R <: AbstractPersistenceReduction, C, FI, IX}
    RM = reduce!(R, boundary(flt, reduced = reduced))
    return PersistentHomology{R,C,FI,IX}(flt, RM, length0)
end

persistenthomology(::Type{PersistentCocycleReduction{R}}, flt::Filtration; length0=false) where {R} =
    PersistentHomology{PersistentCocycleReduction{R}}(flt, Set[], length0)

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

"""Return homology group type: dimension & Betti number"""
eltype(h::PersistentHomology) = Tuple{Int, Int}

function iterate(h::PersistentHomology, p=0)
    p > dim(h.filtration.complex) && return nothing
    βₚ = group(h, p)
    return (p, βₚ), p+1
end
