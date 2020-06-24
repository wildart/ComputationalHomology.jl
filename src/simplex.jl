#=== Simplex ===#
struct Simplex{N, P} <: AbstractSimplex
    vs::NTuple{N, P}
    hash::UInt64
    Simplex{N,P}(vals::NTuple{N, P}) where {N, P} = new(vals, hash(Set(vals)))
    Simplex{0,P}() where {P} = new((), zero(UInt))
end
Simplex(splx::P...) where {P} = Simplex{length(splx),P}(splx)
Simplex(splx::Vector{P}) where {P} = Simplex(splx...)
Simplex(::Type{P}) where {P} = Simplex{0,P}()

# Private methods

show(io::IO, splx::Simplex) = show(io, "σ$(splx.vs)")
eltype(splx::Simplex{N,P}) where {N,P} = P
eltype(::Type{Simplex{N,P}}) where {N,P} = P

# Public methods

dim(splx::Simplex{N,P}) where {N,P} = N-1
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

function faces(σ::Simplex)
    d = dim(σ)
    d == 0 && return Simplex[]
    (Simplex([ifelse(j < i, σ.vs[j], σ.vs[j+1]) for j in 1:d]) for i in 1:(d+1))
end

function boundary(::Type{R}, σ::Simplex) where {R}
    d = dim(σ)
    ch = Chain(d-1, R)
    d == 0 && return ch
    o = one(R)
    for (i,face) in enumerate(faces(σ))
        push!(ch, hash(face)=>( isodd(i) ? o : -o))
    end
    return ch
end

# Public methods for AbstractSimplex

values(splx::Simplex) = collect(splx.vs)


# Misc. methods

function volume(S::AbstractMatrix)
    d, vc = size(S)
    @assert d == vc-1 "Number of vertices in simplex should be dim+1"
    v0 = S[:,1]
    return abs(det(S[:,2:end] .- v0))/prod(1:d)
end
volume(s::Simplex{N,Int}, X::AbstractMatrix) where {N}= volume(X[:, values(s)])

parse(::Type{Simplex{N,P}}, str::AbstractString) where {N,P} =
    Simplex(map(e->parse(P, e), split(str, ' ')))
