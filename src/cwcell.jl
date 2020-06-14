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
