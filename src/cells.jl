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
    hash(cell)

Get a unique identifier of the `cell`
"""
hash(c::AbstractCell) = throw(MethodError(hash,(typeof(c),)))

"""
    union(u, v)

Create a new cell from a combination of vertices of cells `u` and `v`.
"""
union(u::C, v::C) where{C<:AbstractCell} = throw(MethodError(union,(C,C)))

"""Abstract simplex type"""
abstract type AbstractSimplex <: AbstractCell end

"""
    values(cell)

Get an array of `cell` values.
"""
values(c::AbstractSimplex) = throw(MethodError(values,(typeof(c),))) #TODO: consider to move to AbstractSimplex

"""
    vertecies(cell)

Get a collection of `cell` vertecies"""
vertecies(c::AbstractSimplex) = throw(MethodError(vertecies,(typeof(c),)))
