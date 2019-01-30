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
mutable struct Chain{PID} <: AbstractChain
    dim::Int
    coefs::Vector{PID} # coefficients from PID
    elems::Vector{Int} # indexes to elements of E
end
Chain(coefs::Vector{PID}, elems::Vector{Int}) where {PID} = Chain{PID}(0, coefs, elems)
Chain(dim::Int, ::Type{PID}) where {PID} = Chain{PID}(dim, PID[], Int[])
Chain(::Type{PID}) where {PID} = Chain(0,PID)
dim(ch::Chain) = ch.dim
setdim!(ch::Chain, dim::Int) = (ch.dim = dim)

# implement interface
Base.getindex(ch::Chain, i::Integer) = (ch.coefs[i], ch.elems[i])
Base.length(ch::Chain) = length(ch.coefs)
Base.eltype(ch::Chain{PID}) where {PID} = PID

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

function Base.push!(ch::Chain{PID}, coef::PID, idx::Int) where {PID}
    push!(ch.coefs, coef)
    push!(ch.elems, idx)
    return ch
end
Base.push!(ch::Chain{PID}, e::Tuple{PID,Int}) where {PID} = push!(ch, e[1], e[2])
Base.push!(ch::Chain{PID}, e::Pair{PID,Int}) where {PID}  = push!(ch, e[1], e[2])
function Base.append!(a::Chain, b::Chain)
    append!(a.coefs, b.coefs)
    append!(a.elems, b.elems)
    return a
end
+(ch::Chain{PID}, e::Tuple{PID,Int}) where {PID} = push!(ch, e[1], e[2])

+(a::Chain{PID}, b::Chain{PID}) where {PID} = Chain{PID}(dim(a), vcat(a.coefs, b.coefs), vcat(a.elems, b.elems))
-(a::Chain{PID}, b::Chain{PID}) where {PID} = Chain{PID}(dim(a), vcat(a.coefs, -b.coefs), vcat(a.elems, b.elems))

*(ch::Chain{PID}, r::PID) where {PID} = Chain{PID}(dim(ch), ch.coefs*r, ch.elems)
*(r::PID, c::Chain{PID}) where {PID} = c*r

function simplify(ch::Chain)
    PID = eltype(ch)
    tmpChain = Dict{Int,PID}()
    for (c,el) in ch
        if haskey(tmpChain, el)
            tmpChain[el] += c
        else
            tmpChain[el] = c
        end
    end
    cc = Chain(dim(ch), PID)
    for (el,c) in tmpChain
        c != zero(PID) && push!(cc, (c,el))
    end
    return cc
end
