"""Filtration of an abstract complex

We call this sequence of complexes the **filtration** of `f` and
think of it as a construction by adding chunks of simplices at a time `t::FI`.

∅ = K₀ ⊆ K₁ ⊆ . . . ⊆ Kₙ = K.
"""
mutable struct Filtration{C<:AbstractComplex, FI<:AbstractFloat}
    # underlying abstract cell complex
    complex::C
    # total order of simplices as array of (dim, simplex id, filtation value)
    total::Vector{Tuple{Int,UInt,FI}}
    divisions::Number
end

order(flt::Filtration) = flt.total
complex(flt::Filtration) = flt.complex
show(io::IO, flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = print(io, "Filtration($(complex(flt)), $FI)")
eltype(flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = C
valtype(flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = FI
length(flt::Filtration) = isinf(flt.divisions) ? length(unique(e->e[3], order(flt))) : flt.divisions

#
# Constructors
#
Filtration(::Type{C}, ::Type{FI}) where {C <: AbstractComplex, FI} =
    Filtration(C(), Vector{Tuple{Int,UInt,FI}}(), Inf)
Filtration(::Type{C}) where {C <: AbstractComplex} = Filtration(C, Float64)
similar(flt::Filtration{C,FI}) where {C <: AbstractComplex, FI} = Filtration(C, FI)

"""Construct filtration from a cell complex using the order of their appearence in the complex"""
function filtration(cplx::AbstractComplex)
    idx = Vector{Tuple{Int,UInt,Float64}}()
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
    idx = Vector{Tuple{Int,UInt,FI}}()
    for d in 0:dim(cplx)
        for c in cells(cplx, d)
            ci = position(cplx, c)
            push!(idx, (d, hash(c), w[d][ci]))
        end
    end
    sort!(idx, by=x->(x[3], x[1])) # sort by dimension & filtration value
    return Filtration(cplx, idx, divisions)
end

function push!(flt::Filtration{C,FI}, cl::AbstractCell, v::FI; recursive=false) where {C<:AbstractComplex, FI}
    cplx = complex(flt)
    @assert isa(cl, eltype(cplx)) "Complex $(cplx) does not accept $(typeof(cl))"
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
minimum(flt::Filtration) = order(flt)[1][3]
maximum(flt::Filtration) = order(flt)[end][3]

"""
    complex(flt, val)

Return a complex from the filtration `flt` at the filtration value `val`
"""
function complex(flt::Filtration{C,FI}, val::FI) where {C <: AbstractComplex, FI}
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
                push!(bm[i+ridx], revidx[(hash(face), d-1)]+ridx)
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


"""
    similarity(flt::Filtration)

Returns a similarity matrix created from 1-subcomplex of simplicial complex and distance weights of the filtration `flt`.
"""
function similarity(flt::Filtration)
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
function write(io::IO, flt::Filtration)
    cplx = complex(flt)
    for (d, ci, fv) in order(flt)
        for k in values(cplx[ci,d])
            write(io, "$k,")
        end
        write(io, "$fv\n")
    end
end

function read(io::IO, ::Type{Filtration{C,FI}}, ::Type{CL}) where {C <: AbstractComplex, FI, CL<:AbstractCell}
    flt = Filtration(C,FI)
    while !eof(io)
        l = readline(io)
        vals = split(l, ',')
        sval = parse(CL, join(vals[1:end-1], ' '))
        fval = parse(FI, vals[end])
        push!(flt, sval, fval)
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

function writebin(io::IO, flt::Filtration)
    cplx = complex(flt)
    for (d, ci, fv) in order(flt)
        write(io, UInt8(d))
        for k in values(cplx[ci,d])
            write(io, UInt32(k))
        end
        write(io, fv)
    end
end

function readbin(io::IO, ::Type{Filtration{C,FI}}, ::Type{CL}; maxoutdim=1) where {C <: AbstractComplex, FI, CL<:AbstractCell}
    flt = Filtration(C,FI)
    while !eof(io)
        d = read(io, UInt8)
        vals = Int[read(io, UInt32) for i in 1:(d+1)]
        fval = read(io, FI)
        if d <= maxoutdim
            push!(flt, CL(vals), fval)
        end
    end
    return flt
end
