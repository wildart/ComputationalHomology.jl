"""CW complex type"""
mutable struct CWComplex <: AbstractComplex
    cells::Dict{Int,Vector{Cell}}   # cells per dimension
end
show(io::IO, cplx::CWComplex) = print(io, "CWComplex($(size(cplx)))")

(::Type{CWComplex})() = CWComplex(Dict{Int,Vector{Cell}}())
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

"""
    push!(complex, cell)

Attach a new `cell` to to a `complex`.
"""
function push!(cplx::CWComplex, c::Cell; recursive=false)
    cdim = dim(c)
    !haskey(cplx.cells, cdim) && setindex!(cplx.cells, Cell[], cdim)
    push!(cplx.cells[cdim], c)
    return [c]
end
