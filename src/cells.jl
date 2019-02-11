"""Abstract cell type"""
abstract type AbstractCell end

"""Dimension of the cell"""
dim(c::AbstractCell) = throw(MethodError(dim,(typeof(c),)))

"""
    values(cell)

Get an array of `cell` values.
"""
values(c::AbstractCell) = throw(MethodError(values,(typeof(c),)))

"""Cell comparison"""
==(a::AbstractCell, b::AbstractCell) = throw(MethodError(==,(typeof(a), typeof(b))))

"""
    faces(cell)

Get a collection of `cell` faces.
"""
faces(c::AbstractCell) = throw(MethodError(faces,(typeof(c),)))

"""
    vertecies(cell)

Get a collection of `cell` vertecies"""
vertecies(c::AbstractCell) = throw(MethodError(vertecies,(typeof(c),)))

"""Get a unique identifier of the cell"""
hash(c::AbstractCell) = throw(MethodError(hash,(typeof(c),)))

"""
    union(u, v)

Create a new cell from a combination of vertices of cells `u` and `v`.
"""
union(u::C, v::C) where{C<:AbstractCell} = throw(MethodError(union,(C,C)))

abstract type AbstractSimplex <: AbstractCell end

abstract type AbstractHashSimplex <: AbstractSimplex end
