abstract type AbstractChain end

""" k-Chain is a formal linear combination of k-cells.

    Given a set E, we can construct a free R-module M that has E ⊆ M as a basis M, the module M is the module of the formal linear combinations of elements of E, or free module over E, R^{(E)}, s.t. given a finite subset {X_1, ..., X_n} of E, a formal linear combination of X_1, ..., X_n is

    a_1 X_1 + ··· + a_n X_n,

    where the a_i ∈ R.
"""
mutable struct Chain{R} <: AbstractChain
    dim::Int
    coefs::Vector{R}   # coefficients from ring R
    elems::Vector{Int} # indexes to elements of E
end
Chain(coefs::Vector{R}, elems::Vector{Int}) where {R} = Chain{R}(0, coefs, elems)
Chain(dim::Int, ::Type{R}) where {R} = Chain{R}(dim, R[], Int[])
Chain(::Type{R}) where {R} = Chain(0,R)
dim(ch::Chain) = ch.dim
setdim!(ch::Chain, dim::Int) = (ch.dim = dim)

Base.getindex(ch::Chain, i::Integer) = (ch.coefs[i],ch.elems[i])
Base.length(ch::Chain) = length(ch.coefs)
Base.eltype(ch::Chain{R}) where {R} = R
function Base.iterate(ch::Chain, state=1)
    state > length(ch) && return nothing
    return (ch[state], state+1)
end

function Base.show(io::IO, ch::Chain)
    if length(ch.elems) == 0
        print(io, "0")
    else
        for (i,(c,el)) in enumerate(zip(ch.coefs, ch.elems))
            i != 1 && print(io, " + ")
            print(io, "$c[$el]")
        end
    end
end

function Base.push!(ch::Chain{R}, coef::R, idx::Int) where {R}
    push!(ch.coefs, coef)
    push!(ch.elems, idx)
    return ch
end
Base.push!(ch::Chain{R}, e::Tuple{R,Int}) where {R} = push!(ch, e[1], e[2])
Base.push!(ch::Chain{R}, e::Pair{R,Int}) where {R}  = push!(ch, e[1], e[2])
function Base.append!(a::Chain, b::Chain)
    append!(a.coefs, b.coefs)
    append!(a.elems, b.elems)
    return a
end
+(ch::Chain{R}, e::Tuple{R,Int}) where {R} = push!(ch, e[1], e[2])

+(a::Chain{R}, b::Chain{R}) where {R} = Chain{R}(dim(a), vcat(a.coefs, b.coefs), vcat(a.elems, b.elems))
-(a::Chain{R}, b::Chain{R}) where {R} = Chain{R}(dim(a), vcat(a.coefs, -b.coefs), vcat(a.elems, b.elems))

*(ch::Chain{R}, r::R) where {R} = Chain{R}(dim(ch), ch.coefs*r, ch.elems)
*(r::R, c::Chain{R}) where {R} = c*r

function simplify(ch::Chain)
    R = eltype(ch)
    tmpChain = Dict{Int,R}()
    for (c,el) in ch
        if haskey(tmpChain, el)
            tmpChain[el] += c
        else
            tmpChain[el] = c
        end
    end
    cc = Chain(dim(ch), R)
    for (el,c) in tmpChain
        c != zero(R) && push!(cc, (c,el))
    end
    return cc
end
