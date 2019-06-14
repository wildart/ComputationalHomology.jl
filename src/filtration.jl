"""Filtration of an abstract complex

We call this sequence of complexes the **filtration** of `f` and
think of it as a construction by adding chunks of simplices at a time `t::FI`.
∅ = K0 ⊆ K1 ⊆ . . . ⊆ Kn = K.
"""
mutable struct Filtration{C<:AbstractComplex, FI}
    # underlying abstract cell complex
    complex::C
    # total order of simplices as array of (dim, simplex id, filtation value)
    total::Vector{Tuple{Int,Integer,FI}}
    divisions::Number
end

order(flt::Filtration) = flt.total
Base.complex(flt::Filtration) = flt.complex
Base.show(io::IO, flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = print(io, "Filtration($(complex(flt)), $FI)")
Base.valtype(flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = FI
Base.eltype(flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = C
Base.length(flt::Filtration) = isinf(flt.divisions) ? length(unique(e->e[3], order(flt))) : flt.divisions

#
# Constructors
#
Filtration(::Type{C}, ::Type{FI}) where {C <: AbstractComplex, FI} =
    Filtration(C(), Vector{Tuple{Int,Integer,FI}}(), Inf)
Base.similar(flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = Filtration(C, FI)

"""Construct filtration from a cell complex using the order of their appearence in the complex"""
function filtration(cplx::AbstractComplex)
    idx = Vector{Tuple{Int,Integer,Int}}()
    i = 1
    for d in 0:dim(cplx)
        for c in cells(cplx, d)
            push!(idx, (dim(c), hash(c), i))
            i += 1
        end
    end
    return Filtration(cplx, idx, Inf)
end

"""Construct filtration from a cell complex and a complex weight function"""
function filtration(cplx::C, w::Dict{Int,Vector{FI}}; divisions::Number=Inf) where {C<:AbstractComplex, FI}
    idx = Vector{Tuple{Int,Integer,FI}}()
    for d in 0:dim(cplx)
        for c in cells(cplx, d)
            ci = position(cplx, c)
            push!(idx, (d, hash(c), w[d][ci]))
        end
    end
    sort!(idx, by=x->(x[3], x[1])) # sort by dimension & filtration value
    return Filtration(cplx, idx, divisions)
end

function Base.push!(flt::Filtration{C,FI}, cl::AbstractCell, v::FI; recursive=false) where {C<:AbstractComplex, FI}
    cplx = complex(flt)
    @assert isa(cl, celltype(cplx)) "Complex $(cplx) does not accept $(typeof(cl))"
    cls = push!(cplx, cl, recursive=recursive)
    ord = order(flt)
    idx = length(ord) == 0 ? 1 : findlast(e->e[3]<=v, ord)
    for c in sort!(cls, by=s->dim(s))
        if idx == length(ord)
            push!(ord, (dim(c), hash(c), v))
        else
            insert!(ord, idx, (dim(c), hash(c), v))
        end
        idx += 1
    end
    return flt
end

#
# Auxiliary function
#
Base.minimum(flt::Filtration) = order(flt)[1][3]
Base.maximum(flt::Filtration) = order(flt)[end][3]

"""
    complex(flt, val)

Return a complex from the filtration `flt` at the filtration value `val`
"""
function Base.complex(flt::Filtration{C,FI}, val::FI) where {C <: AbstractComplex, FI}
    cplx = complex(flt)
    res = eltype(flt)()
    for (d,si,fv) in order(flt)
        fv > val && break
        push!(res, cplx[si, d])
    end
    return res
end

"""
    boundary(flt::Filtration) -> Vector{<:AbsrtactSet}

Generate a boundary matrix from the filtration `flt` for the persistent homology calculations.
"""
function boundary(flt::Filtration; reduced=false)
    ridx = reduced ? 1 : 0
    # initialize boundary matrix
    cplx = complex(flt)
    sz = sum(size(cplx))+ridx

    # find suitable type for BM storage
    STs = [Int8, Int16, Int32, Int64]
    STi = findfirst(i->i>sz, map(typemax, STs))
    BT = STi == 1 ? BitSet : Set{STs[STi]}

    # create empty BM
    bm = map(i->BT(), 1:sz)

    # fill boundary matrix
    ord = order(flt)
    revidx = Dict((ci, d) => i for (i, (d, ci, fv)) in enumerate(ord))
    for (i, (d, ci, fv)) in enumerate(ord)
        if d > 0
            splx = cplx[ci, d]
            for face in faces(splx)
                fi = cplx[face]
                push!(bm[i+ridx], revidx[(fi, d-1)]+ridx)
            end
        elseif reduced
            push!(bm[i+ridx], 1)
        end
    end
    return bm
end

function sparse(∂::Vector{<:AbstractSet})
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


"""Similarity matrix created from 1-subcomplex of simplicial complex `cplx` and distance weights.
"""
function similarity_matrix(flt::Filtration)
    cplx = complex(flt)
    ord = order(flt)
    C0 = map(hash, cells(cplx, 0))
    N = length(C0)
    @assert N > 0 "Complex should not be empty"
    adj = spzeros(valtype(flt),N,N)
    for c in cells(cplx, 1)
        i, j = map(h->findfirst(isequal(h), C0), vertices(c))
        idx = findfirst(v->v[1]==1 && v[2] == hash(c), ord)
        adj[i, j] = adj[j, i] = ord[idx][3]
    end
    return adj
end


#
# I/O
#
function Base.write(io::IO, flt::Filtration)
    cplx = complex(flt)
    for (d, ci, fv) in order(flt)
        for k in values(cplx[ci,d])
            write(io, "$k,")
        end
        write(io, "$fv\n")
    end
end

function Base.read(io::IO, ::Type{Filtration{C,FI}}) where {C <: AbstractComplex, FI}
    flt = Filtration(C,FI)
    ET = eltype(celltype(complex(flt))())
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
"""Loop through the filtration `flt` producing growing simplicial complexes on every iteration"""
function Base.iterate(flt::Filtration, state=nothing)
    ord = order(flt)
    if state === nothing # calculate initial state
        idx = 1
        fval = ord[idx][3]
        incr = (maximum(flt)-minimum(flt)) / flt.divisions
    else
        idx, fval, incr = state
    end
    idx > length(ord) && return nothing # done
    splxs = Tuple{Int,Integer}[] #simplex dim & index
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

# Filtration simplex iterator
Base.length(splxs::Simplices{<:Filtration}) = length(splxs.itr)
Base.eltype(splxs::Simplices{<:Filtration}) = eltype(splxs.itr)

function Base.iterate(splxs::Simplices{<:Filtration}, state=nothing)
    # call underlying iterator
    res = iterate(splxs.itr, state)
    # final state
    res == nothing && return nothing
    # get complex
    cplx = complex(splxs.itr)
    ss = [cplx[i, d] for (d, i) in res[1][2]]
    return (res[1][1], ss), res[2] # state of filtration iterator
end
