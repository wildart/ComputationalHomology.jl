"""Simplicial Complex Type

Create a simplicial complex with simplices with a value type `P`.
"""
mutable struct SimplicialComplex{S<:AbstractSimplex} <: AbstractComplex{S}
    cells::Dict{Int,Vector{S}}         # cells per dimension
    order::Vector{Dict{UInt64,Int}}  # total order of cells per dimension
end
show(io::IO, cplx::SimplicialComplex) = print(io, "SimplicialComplex($(size(cplx)))")

#
# Constructors
#

SimplicialComplex(::Type{S}) where {S<:AbstractSimplex} =
    SimplicialComplex(Dict{Int,Vector{S}}(), Dict{UInt64,Int}[])
(::Type{SimplicialComplex{S}})() where {S<:AbstractSimplex} = SimplicialComplex(S)

function SimplicialComplex(splxs::S...) where {S<:AbstractSimplex}
    cplx = SimplicialComplex(S)
    added = Set{S}()
    toprocess = Vector{S}()
    for splx in splxs
        push!(toprocess, splx) # add simplex to process its faces
        while(length(toprocess) > 0)
            tmp = pop!(toprocess) # processing simplex faces
            tmp in added && continue # skip if already processed
            addsimplex!(cplx, tmp) # add simples to complex
            push!(added, tmp) # mark as processed
            for f in faces(tmp) # add simplex faces for processing
                push!(toprocess, f)
            end
        end
    end

    return cplx
end

SimplicialComplex(splxs::AbstractVector{<:AbstractSimplex}) = SimplicialComplex(splxs...)

# ---------------
# Private Methods
# ---------------

"""Add a simplex to the complex (as is) and return a complex size, a dimension of simplex and its index in it

This function **doesn't** add missing faces of the simplex to the complex. Use `addsimplices!` function for instead.
"""
function addsimplex!(cplx::SimplicialComplex{S}, splx::S) where {S <: AbstractSimplex}
    d = dim(splx)
    d < 0 && return nothing
    if !haskey(cplx.cells, d)
        cplx.cells[d] = S[]
    end
    if length(cplx.order) <= d
        # println("Missing dimensions: $(length(cplx.order)) need to be $(d+1)")
        for i in length(cplx.order):d
            push!(cplx.order, Dict{UInt64,Int}())
        end
    end
    push!(cplx.cells[d], splx)
    cplx.order[d+1][hash(splx)] = length(cplx.cells[d]) # length(cplx.order)+1
    return splx
end

"""Add a simplex to the complex and all of its faces recursivly and return a complex size, a dimension of simplex and its index in it"""
function addsimplices!(cplx::SimplicialComplex{S}, splx::S) where {S <: AbstractSimplex}
    added = Set{S}() # already cached
    toprocess = S[]
    ret = S[]
    addedidxs = NTuple{3,Integer}[]

    # add simplex to complex
    s = addsimplex!(cplx, splx)
    push!(addedidxs, (sum(size(cplx)), dim(s), hash(s)))
    push!(ret, s)
    for f in faces(splx) # add simplex faces for processing
        push!(toprocess, f)
    end

    # process simplex faces and add missing simplexes
    while(length(toprocess) > 0)
        tmp = pop!(toprocess) # processing simplex faces
        d = dim(tmp)

        tmp ∈ cplx  && continue # skip if already in complex
        tmp ∈ added && continue # skip if already processed

        s = addsimplex!(cplx, tmp) # add simples to the complex
        push!(addedidxs, (sum(size(cplx)), dim(s), hash(s)))
        push!(added, tmp) # mark as processed
        push!(ret, s)

        dim(tmp) == 0 && continue # do not create faces for 0-simplexes
        for f in faces(tmp) # add simplex faces for processing
            push!(toprocess, f)
        end
    end

    return ret
end

#
# AbstractComplex Interface
#
eltype(cplx::SimplicialComplex{S}) where {S<:AbstractSimplex} = S

function cells(cplx::SimplicialComplex{S}) where {S<:AbstractSimplex}
    length(cplx.cells) == 0 && return Vector{S}[] # no cell in complex
    dims = maximum(keys(cplx.cells))
    return [haskey(cplx.cells, d) ? cplx.cells[d] : S[] for d in 0:dims]
end

cells(cplx::SimplicialComplex{S}, d::Int) where {S<:AbstractSimplex} = get(cplx.cells, d,  S[])

function boundary(cplx::SimplicialComplex, idx::IX, d::Int, ::Type{R}) where {R, IX<:Integer}
    ch = Chain(d-1, IX, R)
    d == 0 && return ch

    splx = cplx[idx, d]
    splx === nothing && return ch

    pos = one(R)
    neg = -one(R)
    sgn = true
    for face in faces(splx)
        push!(ch, hash(face)=>(sgn ? pos : neg))
        sgn = !sgn
    end

    return ch
end

function coboundary(cplx::SimplicialComplex, idx::IX, d::Int, ::Type{R}) where {R, IX<:Integer}
    ch = Chain(d+1, IX, R)
    d == 0 && return ch

    # (δϕ)(c) = ϕ(∂c)
    for c in cells(cplx, d+1)
        i = hash(c)
        for (elem, coef) in boundary(cplx, i, d+1, R)
            elem == idx && push!(ch, i=>coef)
        end
    end

    return ch
end

push!(cplx::SimplicialComplex, splx::AbstractSimplex; recursive=false) =
    recursive ? addsimplices!(cplx, splx) : [addsimplex!(cplx, splx)]

function position(cplx::SimplicialComplex, idx::Integer, d::Int)
    length(cplx.order) <= d && return 0
    return get(cplx.order[d+1], idx, 0)
end


#
# Miscellaneous
#

similar(cplx::SimplicialComplex{S}) where {S<:AbstractSimplex} = SimplicialComplex(S)

function read(io::IO, ::Type{SimplicialComplex{S}}) where {S<:AbstractSimplex}
    splxs = S[]
    P = eltype(S())
    while !eof(io)
        l = chomp(readline(io))
        si = map(e->parse(P, e), split(l, ' '))
        push!(splxs, Simplex(si...))
    end
    return SimplicialComplex(splxs)
end

function write(io::IO, cplx::SimplicialComplex)
    zsplxs = Dict{Int,UInt}()
    for s in cells(cplx, 0)
        idx = position(cplx, s)
        zsplxs[first(values(s))] = idx
        write(io, "$idx")
        write(io, 0x0A)
    end
    for d in 1:dim(cplx)
        for s in cplx.cells[d]
            for (j,i) in enumerate(values(s))
                write(io, "$(zsplxs[i])")
                if j>d
                    write(io, 0x0A)
                else
                    write(io, ' ')
                end
            end
        end
    end
    return
end

#
# Complex simplex iterator
#

struct Simplices{T}
    itr::T
    dim::Int
end
simplices(itr::T, dim::Int=-1) where T = Simplices{T}(itr, dim)
show(io::IO, splxs::Simplices{T}) where T =
    print(io, "Simplex Iterator", splxs.dim < 0 ? "" : " (d=$(splxs.dim))", " for $T")

length(splxs::Simplices) = splxs.dim < 0 ? sum(size(splxs.itr)) : size(splxs.itr, splxs.dim)

eltype(iter::Simplices) = eltype(iter.itr)

# State (total id, dim, dim id)
function iterate(iter::Simplices, (tid, d, did)=(0, iter.dim, 1))
    tid >= length(iter) && return nothing
    if d < 0
        d = 0
    end
    if did > size(iter.itr, d)
        d += 1
        did = 1
    end
    return cells(iter.itr, d)[did], (tid+1, d, did+1)
end
