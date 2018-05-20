"""Filtration of an abstract complex

We call this sequence of complexes the **filtration** of `f` and
think of it as a construction by adding chunks of simplices at a time `t::FI`.
∅ = K0 ⊆ K1 ⊆ . . . ⊆ Kn = K.
"""
mutable struct Filtration{C<:AbstractComplex, FI}
    # underlying abstract cell complex
    complex::C
    # total order of simplices as array of (dim, simplex id, filtation value)
    total::Vector{Tuple{Int,Int,FI}}
end
Base.length(flt::Filtration) = length(flt.total)
Base.show(io::IO, flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = print(io, "Filtration($(complex(flt)), $FI)")
Base.valtype(flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = FI

Base.complex(flt::Filtration) = flt.complex
order(flt::Filtration) = flt.total

#
# Constructors
#
Filtration(::Type{C}, ::Type{FI}) where {C <: AbstractComplex, FI} =
    Filtration(C(), Vector{Tuple{Int,Int,FI}}())

"""Construct filtration from a cell complex using the order of their appearence in the complex"""
function filtration(cplx::C) where {C<:AbstractComplex}
    idx = Vector{Tuple{Int,Int,Int}}()
    i = 1
    for d in 0:dim(cplx)
        for c in get(cells(cplx, d))
            push!(idx, (dim(c), c[:index], i))
            i += 1
        end
    end
    return Filtration(cplx, idx)
end

"""Construct filtration from a cell complex and a complex weight function"""
function filtration(cplx::C, w::Dict{Int,Vector{FI}}) where {C<:AbstractComplex, FI}
    idx = Vector{Tuple{Int,Int,FI}}()
    for d in 0:dim(cplx)
        for c in get(cells(cplx, d))
            ci = c[:index]
            push!(idx, (d, ci, w[d][ci]))
        end
    end
    sort!(idx, by=x->(x[3], x[1])) # sort by dimension & filtration value
    return Filtration(cplx, idx)
end

function Base.push!(flt::Filtration{C,FI}, cl::AbstractCell, v::FI; recursive=false) where {C<:AbstractComplex, FI}
    cplx = complex(flt)
    @assert isa(cl, celltype(cplx)) "Complex $(cplx) does not accept $(typeof(cl))"
    cls = push!(cplx, cl, recursive=recursive)
    idx = length(flt.total) == 0 ? 1 : findlast(e->e[3]<=v, flt.total)
    for c in sort!(cls, by=s->dim(s))
        if idx == length(flt.total)
            push!(flt.total, (dim(c), c[:index], v))
        else
            insert!(flt.total, idx, (dim(c), c[:index], v))
        end
        idx += 1
    end
    return flt
end

"""Generate a combined boundary matrix from the filtration `flt` for the persistent homology calculations."""
function boundary_matrix(flt::Filtration; reduced=false)
    ridx = reduced ? 1 : 0
    # initialize boundary matrix
    cplx = complex(flt)
    bm = map(i->IntSet(), 1:sum(size(cplx))+ridx)
    # fill boundary matrix
    for (i, (d, ci, fv)) in enumerate(flt.total)
        if d > 0
            splx = get(cplx[ci, d])
            for face in faces(splx)
                fi = cplx[face, d-1]
                push!(bm[i+ridx], findfirst(e->e[1] == d-1 && e[2] == fi, flt.total)+ridx)
            end
        elseif reduced
            push!(bm[i+ridx], 1)
        end
    end
    return bm
end

function Base.sparse(∂::Vector{IntSet})
    m = length(∂)
    ret = spzeros(Int, m, m)
    for i in 1:m
        bm = ∂[i]
        for (l, j) in enumerate(bm)
            ret[j,i] = j # (-1)^((l-1)%2) # coefs require exact order of faces in provides simplex
        end
    end
    return ret
end

#
# I/O
#

function Base.write(io::IO, flt::Filtration)
    cplx = complex(flt)
    for (d, ci, fv) in flt.total
        for k in get(cplx[ci,d])[:values]
            write(io, "$k,")
        end
        write(io, "$fv\n")
    end
end

function Base.read(io::IO, ::Type{Filtration{C,FI}}) where {C <: AbstractComplex, FI}
    flt = Filtration(C,FI)
    ST = celltype(complex(flt))
    ET = eltype(ST())
    while !eof(io)
        l = readline(io)
        vals = split(l, ',')
        svals = map(v->parse(ET, v), vals[1:end-1])
        fval = parse(FI, vals[end])
        push!(flt, ST(svals), fval)
    end
    return flt
end

"Write a combined boundary matrix to a text file"
function writeboundarymatrix(io::IO, bm::Vector, zeroindex = true)
    for smplx in bm
        if length(smplx) == 0
            write(io, "0")
        else
            write(io, "$(length(smplx)-1)")
        end
        for i in smplx
            write(io, " $(zeroindex ? i-1: i)")
        end
        write(io, 0x0A)
    end
end

#
# Iterator methods
#
function Base.start(flt::Filtration{C, FI}) where {C<:AbstractComplex, FI}
    return (C(), 0)
end

function Base.next(flt::Filtration{C, FI}, state) where {C<:AbstractComplex, FI}
    c = copy(state[1])
    i = state[2]+1
    d, ci, v = flt.total[i]
    push!(c, get(complex(flt)[ci, d]))
    return (v, c), (c, i)
end

function Base.done(flt::Filtration{C, FI}, state) where {C<:AbstractComplex, FI}
     return state[2] == length(flt.total)
end

#
# Filtration simplex iterator
#
Base.length(splxs::Simplices{F}) where {F <: Filtration} = length(unique(e->e[3], order(splxs.itr)))

function Base.start(splxs::Simplices{F}) where {F <: Filtration}
    minval = mapreduce(v->v[3], min, order(splxs.itr))
    return (minval,)
end

function Base.next(splxs::Simplices{F}, state) where {F <: Filtration}
    v = state[1]
    ord = order(splxs.itr)
    cplx = complex(splxs.itr)
    CT = celltype(cplx)
    ss = CT[]
    idx = findfirst(e->e[3] == v, ord)
    if idx == 0
        nextv = Inf
    else
        while length(ord) >= idx && ord[idx][3] == v
            s = cplx[ord[idx][2], ord[idx][1]]
            if !isnull(s)
                push!(ss, get(s))
            end
            idx += 1
        end
        nextv = length(ord) >= idx ? ord[idx][3] : Inf
    end

    return (v, ss), (nextv,)
end

function Base.done(splxs::Simplices{F}, state) where {F <: Filtration}
    return isinf(state[1])
end
