abstract AbstractChain

""" k-Chain is a formal linear combination of k-cells.

    Given a set E, we can construct a free R-module M that has E ⊆ M as a basis M, the module M is the module of the formal linear combinations of elements of E, or free module over E, R^{(E)}, s.t. given a finite subset {X_1, ..., X_n} of E, a formal linear combination of X_1, ..., X_n is

    a_1 X_1 + ··· + a_n X_n,

    where the a_i ∈ R.
"""
type Chain{R} <: AbstractChain
    dim::Int
    coefs::Vector{R}   # coefficients from ring R
    elems::Vector{Int} # indexes to elements of E
end
Chain{R}(coefs::Vector{R}, elems::Vector{Int}) = Chain{R}(0, coefs, elems)
Chain{R}(dim::Int, ::Type{R}) = Chain{R}(dim, R[], Int[])
Chain{R}(::Type{R}) = Chain(0,R)
dim(ch::Chain) = ch.dim
setdim!(ch::Chain, dim::Int) = (ch.dim = dim)

Base.getindex(ch::Chain, i::Integer) = (ch.coefs[i],ch.elems[i])
Base.length(ch::Chain) = length(ch.coefs)
Base.start(ch::Chain) = 1
Base.next(ch::Chain, i) = (ch[i], i+1)
Base.done(ch::Chain, i) = i == length(ch)+1

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

function Base.push!{R}(ch::Chain{R}, coef::R, idx::Int)
    push!(ch.coefs, coef)
    push!(ch.elems, idx)
    return ch
end
Base.push!{R}(ch::Chain{R}, e::Tuple{R,Int}) = push!(ch, e[1], e[2])
Base.push!{R}(ch::Chain{R}, e::Pair{R,Int})  = push!(ch, e[1], e[2])
function Base.append!{R}(a::Chain{R}, b::Chain{R})
    append!(a.coefs, b.coefs)
    append!(a.elems, b.elems)
    return a
end
+{R}(a::Chain{R}, b::Chain{R}) = Chain{R}(dim(a), vcat(a.coefs, b.coefs), vcat(a.elems, b.elems))
-{R}(a::Chain{R}, b::Chain{R}) = Chain{R}(dim(a), vcat(a.coefs, -b.coefs), vcat(a.elems, b.elems))

*{R}(ch::Chain{R}, r::R) = Chain{R}(dim(ch), ch.coefs*r, ch.elems)
*{R}(r::R, c::Chain{R}) = c*r

function simplify{R}(ch::Chain{R})
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
