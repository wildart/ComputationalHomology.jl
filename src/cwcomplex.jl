"""CW complex type"""
mutable struct CWComplex{S<:AbstractCell} <: AbstractComplex{S}
    cells::Dict{Int,Vector{S}}   # cells per dimension
end
show(io::IO, cplx::CWComplex) = print(io, "CWComplex($(size(cplx)))")

(::Type{CWComplex{S}})() where {S<:AbstractCell} = CWComplex{S}(Dict{Int,Vector{S}}())
CWComplex() = CWComplex{Cell{Int}}()
function CWComplex(cells::S...) where {S<:AbstractCell}
    cplx = CWComplex{S}()
    for c in cells
        push!(cplx, c)
    end
    return cplx
end

#
# AbstractComplex Interface
#
eltype(cplx::CWComplex{S}) where {S<:AbstractCell} = S

function cells(cplx::CWComplex{S}) where {S<:AbstractCell}
    length(cplx.cells) == 0 && return Vector{S}[] # no cell in complex
    dims = maximum(keys(cplx.cells))
    return [haskey(cplx.cells, d) ? cplx.cells[d] : S[] for d in 0:dims]
end

cells(cplx::CWComplex{S}, d::Int) where {S<:AbstractCell} = get(cplx.cells, d,  S[])

function boundary(cplx::CWComplex, idx::IX, d::Int, ::Type{PID}) where {PID, IX<:Integer}
end

function coboundary(cplx::CWComplex, idx::IX, d::Int, ::Type{PID}) where {PID, IX<:Integer}
end

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
