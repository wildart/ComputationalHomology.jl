"""CW complex type"""
mutable struct CWComplex <: AbstractComplex
    cells::Dict{Int,Vector{Cell}}   # cells per dimension
    order::Dict{UInt,Tuple{Int,Int}}     # cell index, id => (dim, pos)
end
show(io::IO, cplx::CWComplex) = print(io, "CWComplex($(size(cplx)))")

(::Type{CWComplex})() = CWComplex(Dict{Int,Vector{Cell}}(), Dict{UInt,Tuple{Int,Int}}())
# CWComplex() = CWComplex{Cell{Int}}()
function CWComplex(cells::Cell...)
    cplx = CWComplex()
    for c in cells
        push!(cplx, c)
    end
    return cplx
end

#
# AbstractComplex Interface
#
eltype(cplx::CWComplex) = Cell

function cells(cplx::CWComplex)
    length(cplx.cells) == 0 && return Vector{Cell}[] # no cell in complex
    dims = maximum(keys(cplx.cells))
    return [haskey(cplx.cells, d) ? cplx.cells[d] : Cell[] for d in 0:dims]
end

cells(cplx::CWComplex, d::Int) = get(cplx.cells, d,  Cell[])

function faces(cplx::CWComplex, ci::IX) where {IX<:Integer}
    d, p = cplx.order[ci]
    d == 0 && return IX[]
    return map(hash, faces(cplx.cells[d][p]))
end

function cofaces(cplx::CWComplex, ci::IX) where {IX<:Integer}
    ret = IX[]
    d, p = cplx.order[ci]
    d+=1
    d > dim(cplx) && return ret
    for cid in map(hash, cplx.cells[d])
        ci ∈ faces(cplx, cid) && push!(ret, cid)
    end
    return ret
end

"""
    push!(complex, cell)

Attach a new `cell` to to a `complex`.
"""
function push!(cplx::CWComplex, c::Cell; recursive=false)
    cdim = dim(c)
    !haskey(cplx.cells, cdim) && setindex!(cplx.cells, Cell[], cdim)
    push!(cplx.cells[cdim], c)
    cplx.order[hash(c)] = (cdim, length(cplx.cells[cdim]))
    return [c]
end
