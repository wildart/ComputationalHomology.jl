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
        for c in cells(cplx, d)
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
        for c in cells(cplx, d)
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
    bm = map(i->BitSet(), 1:sum(size(cplx))+ridx)
    # fill boundary matrix
    for (i, (d, ci, fv)) in enumerate(flt.total)
        if d > 0
            splx = cplx[ci, d]
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

function SparseArrays.sparse(∂::Vector{BitSet})
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
        for k in cplx[ci,d][:values]
            write(io, "$k,")
        end
        write(io, "$fv\n")
    end
end

function Base.read(io::IO, ::Type{Filtration{C,FI}}) where {C <: AbstractComplex, FI}
    flt = Filtration(C,FI)
    ET = eltype(celltype(complex(flt)))
    while !eof(io)
        l = readline(io)
        vals = split(l, ',')
        svals = map(v->parse(ET, v), vals[1:end-1])
        fval = parse(FI, vals[end])
        push!(flt, Simplex(svals), fval)
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
            write(io, " $(zeroindex ? i-1 : i)")
        end
        write(io, 0x0A)
    end
end

#
# Iterators
#

# Filtration simplicial complex iterator
Base.length(flt::Filtration) = length(flt.total)
Base.eltype(flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = C

"""Loop through the filtration `flt` producing growing simplicial complexes on every iteration"""
function Base.iterate(flt::Filtration{C, FI}, state=(C(), 0)) where {C <: AbstractComplex, FI}
    state[2] >= length(flt.total) && return nothing # done
    c = copy(state[1]) # form next state copy
    i = state[2]+1
    d, ci, v = flt.total[i]
    push!(c, complex(flt)[ci, d])
    return (v, c), (c, i)
end

# Filtration simplex iterator
Base.length(splxs::Simplices{F}) where F<:Filtration = length(unique(e->e[3], order(splxs.itr)))
Base.eltype(splxs::Simplices{F}) where F<:Filtration = celltype(splxs.itr)

function Base.iterate(splxs::Simplices{F},state=nothing) where F<:Filtration
    # initial state
    state === nothing && return iterate(splxs, mapreduce(v->v[3], min, order(splxs.itr)))
    # final state
    isinf(state) && return nothing
    ord = order(splxs.itr)
    cplx = complex(splxs.itr)
    ss = celltype(cplx)[]
    idx = findfirst(e->e[3] == state, ord)
    nextstate = Inf
    if idx != 0
        while length(ord) >= idx && ord[idx][3] == state
            s = cplx[ord[idx][2], ord[idx][1]]
            if s !== nothing
                push!(ss, s)
            end
            idx += 1
        end
        nextstate = length(ord) >= idx ? ord[idx][3] : Inf
    end

    return (state, ss), nextstate
end

