"""
    EmptyChain

It's a dummy chain class that does not store anything.
"""
struct EmptyChain <: AbstractChain{Nothing,Nothing} end
length(ch::EmptyChain) = 0
dim(ch::EmptyChain) = 0
keys(ch::EmptyChain) = Nothing[]
values(ch::EmptyChain) = Nothing[]
getindex(ch::EmptyChain, i) = nothing


"""
    Chain{IX<:Integer, R}

The chain of integer elements for `R`-module.

This chain implementation uses a dictionary for storing elements and coefficients.
"""
mutable struct Chain{IX<:Integer, R} <: AbstractChain{IX,R}
    dim::Int
    cells::Dict{IX,R}
end

Chain(d::Int, elems::AbstractVector, coefs::AbstractVector) = Chain(d, Dict(zip(elems, coefs)))
Chain(elems::AbstractVector, coefs::AbstractVector) = Chain(0, elems, coefs)
Chain(d::Int, ::Type{IX}, ::Type{R}) where {R, IX<:Integer} = Chain(d, IX[], R[])
Chain(::Type{IX}, ::Type{R}) where {R, IX<:Integer} = Chain(0, IX, R)
Chain(d::Int, ::Type{R}) where {R} = Chain(d, UInt, R)
Chain(::Type{R}) where {R} = Chain(0, UInt, R)

# implement interface
dim(ch::Chain) = ch.dim
length(ch::Chain) = length(ch.cells)
copy(ch::Chain{IX,R}) where {R, IX<:Integer} = Chain{IX,R}(ch.dim, copy(ch.cells))
keys(ch::Chain{IX,R}) where {R, IX<:Integer} = keys(ch.cells)
values(ch::Chain{IX,R}) where {R, IX<:Integer} = values(ch.cells)
getindex(ch::Chain{IX,R}, k::IX) where {R, IX<:Integer} = get(ch.cells, k, zero(R))
setindex!(ch::Chain{IX,R}, c::R, k::IX) where {R, IX<:Integer} = setindex!(ch.cells, c, k)

function push!(ch::Chain{IX,R}, e::Pair{IX,R}) where {R, IX<:Integer}
    (k,v) = e
    if k ∈ ch
        ch.cells[k] += v
    else
        push!(ch.cells, e)
    end
    ch
end

function iterate(ch::Chain{IX,R}, state...) where {R, IX<:Integer}
    y = iterate(ch.cells, state...)
    y === nothing && return nothing
    return (y[1], y[2])
end

function simplify(ch::Chain)
    for (k,v) in ch.cells
        v == 0 && delete!(ch.cells, k)
    end
    return ch
end
