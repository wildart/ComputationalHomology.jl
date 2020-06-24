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
function boundary(::Type{R}, σ::AbstractCell) where {R}
    d = dim(σ)
    ch = Chain(dim(σ)-1, R)
    d == 0 && return ch
    i = one(R)
    sgn = true
    for face in faces(σ)
        push!(ch, hash(face)=>(sgn ? i : -i))
        sgn = !sgn
    end
    return ch
end
boundary(σ::AbstractCell) = boundary(Int, σ)


"""Abstract simplex type"""
abstract type AbstractSimplex <: AbstractCell end

"""
    values(cell)

Get an array of `cell` values.
"""
values(c::AbstractSimplex) = throw(MethodError(values,(typeof(c),)))
