"""Abstract cell type"""
abstract type AbstractCell end
"""Dimension of the cell"""
dim(c::AbstractCell) = throw(MethodError(dim,(typeof(c),)))
"""Get cell properties: `:index` & `:values`"""
Base.getproperty(c::AbstractCell, k::Symbol) = throw(MethodError(getproperty,(typeof(c),)))
"""Set cell properties: `:index` & `:values`"""
Base.setproperty!(c::AbstractCell, v, k::Symbol) = throw(MethodError(setproperty!,(typeof(c),)))
"""Cell comparison"""
==(a::AbstractCell, b::AbstractCell) = throw(MethodError(==,(typeof(a), typeof(b))))
"""Get cell faces"""
faces(c::AbstractCell) = throw(MethodError(faces,(typeof(c),)))

#=== Simples ===#
mutable struct Simplex{P} <: AbstractCell
    idx::Int
    vs::Set{P}
    hash::UInt64
    Simplex{P}(idx, vals) where {P} = new(idx, vals, hash(vals))
    Simplex{P}() where {P} = new(0, Set{P}())
end
Simplex(splx::Vector{P}) where {P} = Simplex{P}(0, Set(splx))
Simplex(splx::P...) where {P} = Simplex(P[splx...])

# Private methods

Base.convert(::Type{Simplex{P}}, v::Vector{P}) where {P} = Simplex{P}(0, Set(v))
Base.hash(splx::Simplex) = splx.hash
Base.show(io::IO, splx::Simplex) = show(io, "σ$(splx.idx)$(collect(splx.vs))")
Base.eltype(splx::Simplex{P}) where {P} = P

# Public methods

dim(splx::Simplex) = length(splx.vs)-1

function Base.getproperty(splx::Simplex, name::Symbol)
    if name == :index
        return splx.idx
    elseif name == :values
        return splx.vs
    else
        return getfield(splx, name)
    end
end

function Base.setproperty!(splx::Simplex, name::Symbol, v)
    if name == :index
        splx.idx = v
    elseif name == :values
        splx.vs = v
    else
        return setfield!(splx, name, v)
    end
end

==(a::Simplex, b::Simplex) = a.hash == b.hash

function faces(splx::Simplex{P}) where {P}
    faces = typeof(splx)[]
    for i in 1:dim(splx)+1
        face = collect(splx.values)
        deleteat!(face,i)
        push!(faces, Simplex(face))
    end
    return faces
end

# Misc. methods

function volume(S::AbstractMatrix)
    d, vc = size(S)
    @assert d == vc-1 "Number of vertexes in simplex should be dim+1"
    v0 = S[:,1]
    return abs(det(S[:,2:end] .- v0))/prod(1:d)
end
volume(s::Simplex{Int}, X::AbstractMatrix) = volume(X[:,collect(s.values)])

#=== Cube ===#

mutable struct Cube{P<:Vector} <: AbstractCell
    origin::P
    extent::P
    function Cube{P}(o::P, x::P) where {P}
        @assert length(x) <= length(o) "Too many extent elements"
        new(o, x)
    end
end
Base.show(io::IO, c::Cube) = show(io, "Cube[$(c.origin) + $(c.extent)]")

# Public methods

dim(c::Cube) = length(findall(!iszero, c.extent))

function Base.getproperty(c::Cube, name::Symbol)
    if name == :index || name == :values
        return (c.origin, c.extent)
    else
        return getfield(c, name)
    end
end

function Base.setproperty!(c::Cube, name::Symbol, v)
    if name == :index || name == :values
        c.origin = v[1]
        c.extent = v[2]
    else
        return setfield!(c, name, v)
    end
end

==(a::Cube, b::Cube) = a.hash == b.hash

# Misc. methods

function volume(c::Cube)
    sizes = similar(c.origin,0)
    for (i,ex) in enumerate(c.extent)
        if ex != 0.
            expoint = copy(c.origin)
            expoint[i] += ex
            l = norm(expoint - c.origin)
            push!(sizes, l)
        end
    end
    return prod(sizes)
end

# Simplex Iterators

struct Simplices{T}
    itr::T
    dim::Int
end
simplices(itr::T, dim::Int=-1) where T = Simplices{T}(itr, dim)
Base.show(io::IO, splxs::Simplices{T}) where T = print(io, "Simplex Iterator", splxs.dim < 0 ? "" : " (d=$(splxs.dim))", " for $T")
