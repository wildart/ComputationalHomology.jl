global cellId = 0

function resetCellId!(id=0)
    global cellId
    cellId = id
end

"""CW cell type"""
struct Cell{R} <: AbstractCell
    d::Int
    boundary::Vector{Cell{R}}
    boundaryMap::Dict{UInt64,R}
    hash::UInt64
    function Cell{R}(d::Int, cells::Vector{Cell{R}}, boundaryMap::Dict{UInt64,R}) where {R}
        global cellId
        @assert length(cells) == length(boundaryMap) "Collections of cells and boundary coefficients must be of the same size"
        cellId += 1
        new(d, cells, boundaryMap, hash(cellId))
    end
end
function Cell(d::Int, cells::Vector{Cell{R}}, boundary::Vector{R}) where {R}
    boundaryMap = Dict{UInt64,R}( hash(c)=>b for (c,b) in zip(cells, boundary) )
    Cell{R}(d, cells, boundaryMap)
end
function Cell(d::Int, cells::Vector{Cell{R}}) where {R}
    o = one(R)
    boundaryMap = Dict{UInt64,R}()
    for (i,c) in enumerate(cells)
        @assert dim(c) == d-1 "Incorrect dimension of a boundary cell $c"
        boundaryMap[hash(c)] = isodd(i) ? o : -o
    end
    Cell{R}(d, cells, boundaryMap)
end
Cell(d::Int, cells::Cell{R}...) where {R} = Cell(d, [cells...])
Cell(::Type{R}, d::Int) where {R} = Cell(d, Cell{R}[])
Cell() = Cell(Int, 0)

# Private methods

function show(io::IO, c::Cell)
    vstr = ""
    print(io, "C$(c.d)")
    if length(c.boundary) > 0
        vstr = foldl((v,u)->v*"$(u){$(c.boundaryMap[hash(u)])}, ", c.boundary, init="")[1:end-2]
        print(io, "[$vstr]")
    end
end
eltype(c::Cell{R}) where {R} = R
eltype(::Type{Cell{R}}) where {R} = R

function push!(c1::Cell{R}, c2::Cell{R}, coeff::R) where {R}
    @assert dim(c1)-1 == dim(c2) "Incorrect dimension of a boundary cell $c2"
    h = hash(c2)
    c1.boundaryMap[h] = coeff
    push!(c1.boundary, c2)
end

function push!(c1::Cell{R}, c2::Cell{R}) where {R}
    v = isodd(length(c1.boundaryMap)) ? -one(R) : one(R)
    push!(c1, c2, v)
    c1
end


# Public methods

dim(c::Cell) = c.d
hash(c::Cell) = c.hash
==(a::Cell, b::Cell) = a.hash == b.hash

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

faces(c::Cell) = c.boundary

boundary(c::Cell) = Chain(dim(c), copy(c.boundaryMap))
