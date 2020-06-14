#=== Simplex ===#
struct Simplex{P} <: AbstractSimplex
    vs::Set{P}
    hash::UInt64
    Simplex{P}(vals) where {P} = new(vals, hash(vals))
    Simplex{P}() where {P} = new(Set{P}())
end
Simplex(splx::Vector{P}) where {P} = Simplex{P}(Set(splx))
Simplex(splx::P...) where {P} = Simplex(P[splx...])

# Private methods

Base.convert(::Type{Simplex{P}}, v::Vector{P}) where {P} = Simplex{P}(0, Set(v))
show(io::IO, splx::Simplex) = show(io, "σ$(collect(splx.vs))")
eltype(splx::Simplex{P}) where {P} = P

# Public methods

dim(splx::Simplex) = length(splx.vs)-1

hash(splx::Simplex) = splx.hash

==(a::Simplex, b::Simplex) = a.hash == b.hash

function vertices(splx::Simplex)
    vtxs = typeof(hash(splx))[]
    for v in splx.vs
        push!(vtxs, hash(Simplex(v)))
    end
    return vtxs
end

union(u::Simplex, v::Simplex) = Simplex(collect(u.vs ∪ v.vs))

# Public methods for Simplex

function faces(splx::Simplex)
    faces = typeof(splx)[]
    for i in 1:dim(splx)+1
        face = collect(splx.vs)
        deleteat!(face,i)
        push!(faces, Simplex(face))
    end
    return faces
end

values(splx::Simplex) = collect(splx.vs)

# Misc. methods

function volume(S::AbstractMatrix)
    d, vc = size(S)
    @assert d == vc-1 "Number of vertices in simplex should be dim+1"
    v0 = S[:,1]
    return abs(det(S[:,2:end] .- v0))/prod(1:d)
end
volume(s::Simplex{Int}, X::AbstractMatrix) = volume(X[:,collect(s.vs)])
