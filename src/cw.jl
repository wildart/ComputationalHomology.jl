#=== CW cell ===#
struct Cell <: AbstractCell
    d::Int
    boundary::Vector{Cell}
    hash::UInt64
    Cell(d::Int, cells::Vector{Cell}) = new(d, cells, hash(cells))
end
Cell() = Cell(0, Cell[])
Cell(d::Int, splx::Cell...)= Cell(d, [splx...])

# Private methods

function Base.show(io::IO, c::Cell)
    vstr = ""
    print(io, "C$(c.d)")
    if length(c.boundary) > 0
        vstr = foldl((v,u)->v*"$(u), ", c.boundary, init="")[1:end-2]
        print(io, "[$vstr]")
    end
end
Base.eltype(c::Cell) = Cell

# Public methods

dim(c::Cell) = c.d

==(a::Cell, b::Cell) = a.hash == b.hash

hash(c::Cell) = c.hash

faces(c::Cell) = c.boundary

function vertecies(c::Cell)
    ret = Cell[]
    return ret
end

function union(u::Cell, v::Cell)
    dim(u) != dim(v) && throw(DimensionMismatch("Cell dimensions must match for union"))
    return Cell(dim(v)+1, u, v)
end
