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
    iterate(flt::Filtration{C,FI}) -> Union{Nothing, Tuple{FI,Vector{Tuple{Int,UInt}}}}

Advance the filtration iterator `flt` to obtain the next filtration value of type `FI` with
a corresponding collection of cell parameters: dimension & identifier.
"""
function iterate(flt::Filtration{C,FI}, state=nothing) where {C<:AbstractComplex, FI<:AbstractFloat}
    ord = order(flt)
    if state === nothing # calculate initial state
        idx = 1
        fval = ord[idx][3]
        incr = (maximum(flt)-minimum(flt)) / flt.divisions
    else
        idx, fval, incr = state
    end
    idx > length(ord) && return nothing # done
    splxs = Tuple{Int,UInt}[] #simplex dim & index
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
struct PersistentHomology{R<:AbstractPersistenceReduction, T} <: AbstractHomology
    diagram::Dict{Int, PersistenceDiagram{T}}
end
function persistenthomology(::Type{R}, flt::Filtration{C,FI};
                            kwargs...) where {R<:AbstractPersistenceReduction, C<:AbstractComplex, FI}
    dgm = diagram(R, flt; kwargs...)
    return PersistentHomology{R,FI}(dgm)
end
persistenthomology(flt::Filtration) = persistenthomology(TwistReduction, flt)

diagram(ph::PersistentHomology) = ph.diagram
show(io::IO, h::PersistentHomology{R}) where {R <: AbstractPersistenceReduction} =
    print(io, "Persistent Homology with $R")

#
# Iterator methods
#
values(ph::PersistentHomology) = Iterators.flatten( [birth.(dgm); death.(dgm)] for (d,dgm) in ph.diagram) |> unique |> sort!
length(ph::PersistentHomology) = length(values(ph))

"""Return persistent homology group type: filtraton value & Betti numbers"""
eltype(ph::PersistentHomology{R,T}) where {R, T} = Tuple{T, NTuple}

function iterate(ph::PersistentHomology,
                 p=(1, sort!(collect(keys(ph.diagram))), values(ph)) )
    fidx, dims, fvals = p
    fidx > length(fvals) && return nothing
    maxd = length(dims)
    fv = fvals[fidx]
    βₚ = ntuple(maxd) do i
        if isinf(fv)
            count(v->isinf(death(v)), ph.diagram[i-1])
        else
            count(v->fv ∈ v, ph.diagram[i-1])
        end
    end
    return (fv, βₚ), (fidx+1, dims, fvals)
end

betti(g::PersistentHomology) = map(last, g)

generators(ph::PersistentHomology{PersistentCocycleReduction{R}, T}) where {R, T} =
    Dict( d=>getfield.(cyc, :g) for (d, cyc) in ph.diagram)
