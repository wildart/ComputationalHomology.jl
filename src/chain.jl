abstract type AbstractChain end

# Chain interface
Base.getindex(ch::AbstractChain, i::Integer) = throw(MethodError(getindex, (typeof(ch),Integer)))
Base.length(ch::AbstractChain) = throw(MethodError(length, (typeof(ch),)))
Base.eltype(ch::AbstractChain) = throw(MethodError(eltype, (typeof(ch),)))
function Base.iterate(ch::AbstractChain, state=1)
    state > length(ch) && return nothing
    return (ch[state], state+1)
end

"""Empty Chain"""
struct EmptyChain <: AbstractChain end
Base.length(ch::EmptyChain) = 0

""" k-Chain is a formal linear combination of k-cells.

    Given a set E, we can construct a free R-module M that has E ⊆ M as a basis M, the module M is the module of the formal linear combinations of elements of E, or free module over E, R^{(E)}, s.t. given a finite subset {X_1, ..., X_n} of E, a formal linear combination of X_1, ..., X_n is

    a_1 X_1 + ··· + a_n X_n,

    where the a_i ∈ R.
"""
mutable struct Chain{PID, IX<:Integer} <: AbstractChain
    dim::Int
    coefs::Vector{PID} # coefficients from PID
    elems::Vector{IX}  # indexes to elements of E
end
Chain(coefs::Vector{PID}, elems::Vector{IX}) where {PID, IX<:Integer} = Chain{PID, IX}(0, coefs, elems)
Chain(dim::Int, ::Type{PID}, ::Type{IX}) where {PID, IX<:Integer} = Chain{PID, IX}(dim, PID[], IX[])
Chain(dim::Int, ::Type{PID}) where {PID} = Chain(dim, PID, Int)
Chain(::Type{PID}) where {PID} = Chain(0,PID)
dim(ch::Chain) = ch.dim

# implement interface
Base.getindex(ch::Chain, i::Integer) = (ch.coefs[i], ch.elems[i])
Base.length(ch::Chain) = length(ch.coefs)
Base.eltype(ch::Chain{PID, IX}) where {PID, IX<:Integer} = (PID, IX)
Base.iszero(ch::Chain) = length(ch.elems) == 0

function Base.show(io::IO, ch::Chain)
    if iszero(ch)
        print(io, "0")
    else
        for (i,(c,el)) in enumerate(zip(ch.coefs, ch.elems))
            i != 1 && print(io, " + ")
            print(io, "$c[$el]")
        end
    end
end

function Base.push!(ch::Chain{PID, IX}, coef::PID, idx::IX) where {PID, IX<:Integer}
    push!(ch.coefs, coef)
    push!(ch.elems, idx)
    return ch
end
Base.push!(ch::Chain{PID, IX}, e::Tuple{PID,IX}) where {PID, IX<:Integer} = push!(ch, e[1], e[2])
Base.push!(ch::Chain{PID, IX}, e::Pair{PID,IX}) where {PID, IX<:Integer}  = push!(ch, e[1], e[2])
function Base.append!(a::Chain, b::Chain)
    append!(a.coefs, b.coefs)
    append!(a.elems, b.elems)
    return a
end
+(ch::Chain{PID}, e::Tuple{PID,Integer}) where {PID} = push!(ch, e[1], e[2])

+(a::Chain{PID, IX}, b::Chain{PID, IX}) where {PID, IX<:Integer} = Chain{PID, IX}(dim(a), vcat(a.coefs, b.coefs), vcat(a.elems, b.elems))
-(a::Chain{PID, IX}, b::Chain{PID, IX}) where {PID, IX<:Integer} = Chain{PID, IX}(dim(a), vcat(a.coefs, -b.coefs), vcat(a.elems, b.elems))

*(ch::Chain{PID, IX}, r::PID) where {PID, IX<:Integer} = Chain{PID, IX}(dim(ch), ch.coefs*r, ch.elems)
*(r::PID, c::Chain{PID}) where {PID} = c*r

function simplify(ch::Chain)
    PID, IX = eltype(ch)
    tmpChain = Dict{Integer,PID}()
    for (c,el) in ch
        if haskey(tmpChain, el)
            tmpChain[el] += c
        else
            tmpChain[el] = c
        end
    end
    cc = Chain(dim(ch), PID, IX)
    for (el,c) in tmpChain
        c != zero(PID) && push!(cc, (c,el))
    end
    return cc
end
