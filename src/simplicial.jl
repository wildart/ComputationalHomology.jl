#=== Simplex ===#
struct Simplex{P} <: AbstractSimplex
    vs::Set{P}
    hash::UInt64
    Simplex{P}(vals) where {P} = new(vals, hash(vals))
    Simplex{P}() where {P} = new(Set{P}())
end
Simplex(splx::Vector{P}) where {P} = Simplex{P}(Set(splx))
Simplex(splx::P...) where {P} = Simplex(P[splx...])

# Private methods

Base.convert(::Type{Simplex{P}}, v::Vector{P}) where {P} = Simplex{P}(0, Set(v))
Base.hash(splx::Simplex) = splx.hash
Base.show(io::IO, splx::Simplex) = show(io, "σ$(collect(splx.vs))")
Base.eltype(splx::Simplex{P}) where {P} = P

# Public methods

dim(splx::Simplex) = length(splx.vs)-1

function Base.getproperty(splx::Simplex, name::Symbol)
    if name == :index
        return hash(splx)
    elseif name == :values
        return splx.vs
    else
        return getfield(splx, name)
    end
end

==(a::Simplex, b::Simplex) = a.hash == b.hash

function faces(splx::Simplex)
    faces = typeof(splx)[]
    for i in 1:dim(splx)+1
        face = collect(splx.values)
        deleteat!(face,i)
        push!(faces, Simplex(face))
    end
    return faces
end

function vertecies(splx::Simplex)
    vertecies = typeof(hash(splx))[]
    for v in splx.values
        push!(vertecies, hash(Simplex(v)))
    end
    return vertecies
end

union(u::Simplex, v::Simplex) = Simplex(collect(u.values ∪ v.values))

# Misc. methods

function volume(S::AbstractMatrix)
    d, vc = size(S)
    @assert d == vc-1 "Number of vertexes in simplex should be dim+1"
    v0 = S[:,1]
    return abs(det(S[:,2:end] .- v0))/prod(1:d)
end
volume(s::Simplex{Int}, X::AbstractMatrix) = volume(X[:,collect(s.values)])


# iterator

struct Simplices{T}
    itr::T
    dim::Int
end
simplices(itr::T, dim::Int=-1) where T = Simplices{T}(itr, dim)
Base.show(io::IO, splxs::Simplices{T}) where T = print(io, "Simplex Iterator", splxs.dim < 0 ? "" : " (d=$(splxs.dim))", " for $T")


# Simplicial Complex

"""Simplicial Complex Type

Create a simplicial complex with simplices with a value type `P`.
"""
mutable struct SimplicialComplex{S<:AbstractSimplex} <: AbstractComplex
    cells::Dict{Int,Vector{S}}   # cells per dimension
    order::Dict{UInt64,Int}      # total order of cells per dimension
end
Base.show(io::IO, cplx::SimplicialComplex) = print(io, "SimplicialComplex($(size(cplx)))")

#
# Constructors
#

SimplicialComplex(::Type{S}) where {S<:AbstractSimplex} = SimplicialComplex(Dict{Int,Vector{S}}(), Dict{UInt64,Int}())
(::Type{SimplicialComplex{S}})() where {S<:AbstractSimplex} = SimplicialComplex(S)


function SimplicialComplex(splxs::Simplex...)
    S = eltype(splxs)
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

SimplicialComplex(splxs::Vector{<:AbstractSimplex}) = SimplicialComplex(splxs...)

Base.similar(cplx::SimplicialComplex) = SimplicialComplex(eltype(celltype(cplx)))

# ---------------
# Private Methods
# ---------------

"""Add a simplex to the complex (as is) and return a complex size, a dimension of simplex and its index in it

This function **doesn't** add missing faces of the simplex to the complex. Use `addsimplices!` function for instead.
"""
function addsimplex!(cplx::SimplicialComplex, splx::AbstractSimplex)
    d = dim(splx)
    d < 0 && return nothing
    !haskey(cplx.cells, d) && setindex!(cplx.cells, celltype(cplx)[], d)
    push!(cplx.cells[d], splx)
    cplx.order[hash(splx)] = length(cplx.cells[d]) # length(cplx.order)+1
    return splx
end

"""Add a simplex to the complex and all of its faces recursivly and return a complex size, a dimension of simplex and its index in it"""
function addsimplices!(cplx::SimplicialComplex, splx::AbstractSimplex)
    CT = celltype(cplx)
    # cache already
    added = Set{CT}()

    toprocess = CT[]
    ret = CT[]
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

        cplx[tmp] !== nothing && continue # skip if already in complex
        tmp in added && continue # skip if already processed

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
# Public Interface
#
celltype(cplx::SimplicialComplex) = eltype(valtype(cplx.cells))

function cells(cplx::SimplicialComplex)
    CCT = valtype(cplx.cells)
    length(cplx.cells) == 0 && return CCT[] # no cell in complex
    dims = maximum(keys(cplx.cells))
    return CCT[haskey(cplx.cells, d) ? cplx.cells[d] : CCT() for d in 0:dims]
end

cells(cplx::SimplicialComplex, d::Int) = get(cplx.cells, d,  nothing)

function boundary(cplx::SimplicialComplex, idx::IX, d::Int, ::Type{PID}) where {PID, IX<:Integer}
    ch = Chain(d-1, PID, IX)
    d == 0 && return ch

    pos = one(PID)
    neg = -one(PID)

    splx = cplx[idx, d]
    splx === nothing && return ch

    sgn = true
    for face in faces(splx)
        fidx = cplx[face]
        push!(ch, (sgn ? pos : neg), fidx)
        sgn = !sgn
    end

    return ch
end

function coboundary(cplx::SimplicialComplex, idx::IX, d::Int, ::Type{PID}) where {PID, IX<:Integer}
    ch = Chain(d+1, PID, IX)
    d == 0 && return ch

    # (δϕ)(c) = ϕ(∂c)
    for c in cells(cplx, d+1)
        i = hash(c)
        for (coef, elem) in boundary(cplx, i, d+1, PID)
            elem == idx && push!(ch, coef, i)
        end
    end

    return ch
end

Base.push!(cplx::SimplicialComplex, splx::AbstractSimplex; recursive=false) =
    recursive ? addsimplices!(cplx, splx) : [addsimplex!(cplx, splx)]

Base.position(cplx::SimplicialComplex, idx::Integer, d::Int) where {C<:AbstractCell} = get(cplx.order, idx, nothing)

#
# Miscellaneous
#

function Base.read(io::IO, ::Type{SimplicialComplex{S}}) where {S<:AbstractSimplex}
    splxs = S[]
    P = eltype(S())
    while !eof(io)
        l = chomp(readline(io))
        si = map(e->parse(P, e), split(l, ' '))
        push!(splxs, Simplex(si...))
    end
    return SimplicialComplex(splxs)
end

function Base.write(io::IO, cplx::SimplicialComplex)
    zsplxs = Dict{Int,UInt}()
    for s in cplx.cells[0]
        idx = cplx.order[hash(s)]
        zsplxs[first(s.values)] = idx
        write(io, "$idx")
        write(io, 0x0A)
    end
    for d in 1:dim(cplx)
        for s in cplx.cells[d]
            for (j,i) in enumerate(s.values)
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
Base.length(splxs::Simplices{C}) where C<:AbstractComplex =
    splxs.dim < 0 ? sum(size(splxs.itr)) : size(splxs.itr, splxs.dim)

Base.eltype(iter::Simplices{C}) where C<:AbstractComplex = celltype(iter.itr)

# State (total id, dim, dim id)
function Base.iterate(iter::Simplices{C}, (tid, d, did)=(0, iter.dim, 1)) where C<:AbstractComplex
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
