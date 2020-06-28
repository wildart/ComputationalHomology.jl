"""Abstract cell type"""
abstract type AbstractCell end

"""Dimension of the cell"""
dim(c::AbstractCell) = throw(MethodError(dim,(typeof(c),)))

"""Cell comparison"""
==(a::AbstractCell, b::AbstractCell) = throw(MethodError(==,(typeof(a), typeof(b))))

"""
    faces(cell)

Get a collection of `cell` faces.
"""
faces(c::AbstractCell) = throw(MethodError(faces,(typeof(c),)))

"""
    hash(cell) -> UInt

Return a unique identifier of the `cell` of `UInt` type
"""
hash(c::AbstractCell) = throw(MethodError(hash,(typeof(c),)))

"""
    union(u, v)

Create a new cell from a combination of vertices of cells `u` and `v`.
"""
union(u::C, v::C) where {C<:AbstractCell} = throw(MethodError(union,(C,C)))

"""
    vertices(cell)

Get a collection of `cell` vertices.
"""
vertices(c::AbstractCell) = throw(MethodError(vertices,(typeof(c),)))

"""
    boundary(cell)

Calculate boundary chain of the `cell`.
"""
boundary(::Type{R}, c::AbstractCell) where {R} = throw(MethodError(boundary,(typeof(c),)))
boundary(σ::AbstractCell) = boundary(Int, σ)


(+)(ch::AbstractChain{IX,R}, e::Tuple{C,R}) where {C<:AbstractCell,IX,R} = ch + (hash(e[1]), e[2])


"""Abstract simplex type"""
abstract type AbstractSimplex <: AbstractCell end

"""
    values(cell)

Get an array of `cell` values.
"""
values(c::AbstractSimplex) = throw(MethodError(values,(typeof(c),)))
