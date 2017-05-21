import Base: ==

"""Abstract cell type"""
abstract AbstractCell
"""Dimension of the cell"""
dim(c::AbstractCell) = throw(MethodError(dim,(typeof(c),)))
"""Get cell properties"""
Base.getindex(c::AbstractCell, k::Symbol) = throw(MethodError(getindex,(typeof(c),)))
"""Set cell properties"""
Base.setindex!(c::AbstractCell, v, k::Symbol) = throw(MethodError(setindex!,(typeof(c),)))
"""Cell comparison"""
==(a::AbstractCell, b::AbstractCell) = throw(MethodError(==,(typeof(a), typeof(b))))
"""Get cell faces"""
faces(c::AbstractCell) = throw(MethodError(faces,(typeof(c),)))

#=== Simples ===#
type Simplex{P} <: AbstractCell
    idx::Int
    vs::Vector{P}
end
Simplex{P}(splx::P...) = Simplex{P}(0, sort!([i for i in splx]))
Simplex{P}(splx::Vector{P}) = Simplex{P}(0, sort!(splx))

# Private methods

Base.convert{P}(::Type{Simplex{P}}, v::Vector{P}) = Simplex{P}(0, sort!(v))
Base.hash(splx::Simplex) = hash(splx.vs)
Base.show(io::IO, splx::Simplex) = show(io, "Σ($(splx.vs))[$(splx.idx)]")
Base.eltype{P}(splx::Simplex{P}) = P

# Public methods

dim(splx::Simplex) = length(splx.vs)-1

function Base.getindex(splx::Simplex, k::Symbol)
    if k == :index
        return splx.idx
    elseif k == :values
        return splx.vs
    else
        throw(KeyError(k))
    end
end

function Base.setindex!(splx::Simplex, v, k::Symbol)
    if k == :index
        splx.idx = v
    elseif k == :values
        splx.vs = v
    else
        throw(KeyError(k))
    end
end

==(a::Simplex, b::Simplex) = a[:values] == b[:values]

function faces{P}(splx::Simplex{P})
    faces = typeof(splx)[]
    for i in 1:dim(splx)+1
        face = copy(splx[:values])
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
volume(s::Simplex{Int}, X::AbstractMatrix) = volume(X[:,s[:values]])

#=== Cube ===#

type Cube{P<:Vector} <: AbstractCell
    origin::P
    extent::P
    function Cube(o::P, x::P)
        @assert length(x) <= length(o) "Too many extent elements"
        new(o, x)
    end
end
Base.show(io::IO, c::Cube) = show(io, "Cube[$(c.origin) + $(c.extent)]")

# Public methods

dim(c::Cube) = length(findn(c.extent))

function Base.getindex(c::Cube, k::Symbol)
    if k == :index || k == :values
        return (c.origin, c.extent)
    else
        throw(KeyError(k))
    end
end

function Base.setindex!(c::Cube, v, k::Symbol)
    if k == :index || k == :values
        c.origin = v[1]
        c.extent = v[2]
    else
        throw(KeyError(k))
    end
end

==(a::Cube, b::Cube) = a[:values] == b[:values]

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
