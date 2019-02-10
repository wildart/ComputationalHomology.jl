"""Abstract cell type"""
abstract type AbstractCell end

"""Dimension of the cell"""
dim(c::AbstractCell) = throw(MethodError(dim,(typeof(c),)))

"""Get cell properties: `index` & `values`"""
Base.getproperty(c::AbstractCell, k::Symbol) = throw(MethodError(getproperty,(typeof(c),)))

"""Cell comparison"""
==(a::AbstractCell, b::AbstractCell) = throw(MethodError(==,(typeof(a), typeof(b))))

"""Get cell faces"""
faces(c::AbstractCell) = throw(MethodError(faces,(typeof(c),)))

"""Get cell faces"""
vertecies(c::AbstractCell) = throw(MethodError(vertecies,(typeof(c),)))

"""Get a unique identifier of the cell"""
Base.hash(c::AbstractCell) = throw(MethodError(hash,(typeof(c),)))

Base.union(u::C, v::C) where{C<:AbstractCell} = throw(MethodError(union,(C,C)))

abstract type AbstractSimplex <: AbstractCell end

abstract type AbstractHashSimplex <: AbstractSimplex end
