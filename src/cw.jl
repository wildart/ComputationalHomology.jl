global cellId = 0

function resetCellId!(id=0)
    global cellId
    cellId = id
end

"""CW cell type"""
struct Cell{PID} <: AbstractCell
    d::Int
    boundary::Vector{Cell}
    boundarymap::Dict{UInt64,PID}
    hash::UInt64
    function Cell{PID}(d::Int, cells::Vector{Cell{PID}}) where {PID}
        global cellId
        o = one(PID)
        boundaryMap = Dict{UInt64,PID}()
        for (i,c) in enumerate(cells)
            @assert dim(c) == d-1 "Incorrect dimension of a boundary cell $c"
            h = hash(c)
            boundaryMap[h] = iseven(i) ? -o : o
        end
        cellId += 1
        new(d, cells, boundaryMap, hash(cellId))
    end
end
Cell() = Cell{Int}(0, Cell{Int}[])
Cell(d::Int, cells::Cell{PID}...) where {PID} = Cell{PID}(d, [cells...])
function Cell(d::Int, cells::Vector{Cell{PID}}, boundary::Vector{PID}) where {PID}
    @assert length(cells) == length(boundary) "Collections of cells and boundary coefficients must be of the same size"
    nc = Cell{PID}(d, [cells...])
    for (i,c) in enumerate(nc.boundary)
        nc.boundarymap[hash(c)] = boundary[i]
    end
    return nc
end

# Private methods

function show(io::IO, c::Cell)
    vstr = ""
    print(io, "C$(c.d)")
    if length(c.boundary) > 0
        vstr = foldl((v,u)->v*"$(u){$(c.boundarymap[hash(u)])}, ", c.boundary, init="")[1:end-2]
        print(io, "[$vstr]")
    end
end
eltype(c::Cell{PID}) where {PID} = Cell{PID}

# Public methods

dim(c::Cell) = c.d

==(a::Cell, b::Cell) = a.hash == b.hash

hash(c::Cell) = c.hash

faces(c::Cell) = c.boundary

function vertices(c::Cell)
    vs = Set{Cell}()
    process = Set{Cell}(c.boundary)
    dim(c) == 0 && return collect(vs)
    while length(process) > 0
        pc = pop!(process)
        if dim(pc) == 0
            push!(vs, pc)
        else
            for v in pc.boundary
                push!(process, v)
            end
        end
    end
    return collect(vs)
end

function union(u::Cell, v::Cell)
    dim(u) != dim(v) && throw(DimensionMismatch("Cell dimensions must match for union"))
    return Cell(dim(v)+1, u, v)
end

# Private methods

boundary(c::Cell) = Chain(dim(c), collect(values(c.boundarymap)), collect(keys(c.boundarymap)))

"""Cell complex type"""
mutable struct CWComplex{S<:AbstractCell} <: AbstractComplex
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
celltype(cplx::CWComplex{S}) where {S} = S

function cells(cplx::CWComplex)
    CCT = valtype(cplx.cells)
    length(cplx.cells) == 0 && return CCT[] # no cell in complex
    dims = maximum(keys(cplx.cells))
    return CCT[haskey(cplx.cells, d) ? cplx.cells[d] : CCT() for d in 0:dims]
end

cells(cplx::CWComplex, d::Int) = get(cplx.cells, d,  nothing)

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
