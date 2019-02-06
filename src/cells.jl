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
