"""
    AbstractChain{C, R}

A chain is a formal linear combination of `C`-type elements in `R`-module.

Given a `C`-type element set ``E``, we can construct a free `R`-module ``M`` that
has ``E ⊆ M`` as a basis ``M``, the module ``M`` is the module of the formal linear combinations
of elements of ``E``, or free module over ``E``, ``R^{(E)}``, s.t. given a finite subset
``{X_1, ..., X_n}`` of ``E``, a *formal linear combination* of ``X_1, ..., X_n`` is

``a_1 X_1 + ··· + a_n X_n``

where the ``a_i`` ∈ `R`.
"""
abstract type AbstractChain{C,R} end

# Chain interface
"""
    dim(ch::AbstractChain)

Returns the dimension of the chain `ch`.
"""
dim(ch::AbstractChain) = throw(MethodError(dim, (typeof(ch),)))

"""
    length(ch::AbstractChain)

Returns the length of the chain `ch`.
"""
length(ch::AbstractChain) = throw(MethodError(length, (typeof(ch),)))

"""
    keys(ch::AbstractChain)

Returns elements of the chain `ch`.
"""
keys(ch::AbstractChain) = throw(MethodError(keys, (typeof(ch),)))

"""
    values(ch::AbstractChain)

Returns coefficients of the chain `ch`.
"""
values(ch::AbstractChain) = throw(MethodError(values, (typeof(ch),)))

"""
    copy(ch::AbstractChain)

Returns the copy of the chain `ch`.
"""
copy(ch::AbstractChain) = throw(MethodError(copy, (typeof(ch),)))

"""
    setindex!(ch::AbstractChain{C,R}, v::R, e::C)

Sets an element `e` coefficient `v` in the chain `ch`.
"""
setindex!(ch::AbstractChain{C,R}, v::R, e::C) where {C,R} = throw(MethodError(setindex!, (typeof(ch), R, C)))

"""
    getindex(ch::AbstractChain{C,R}, e::C)

Returns a coefficient of the element `e` of the chain `ch`.
"""
getindex(ch::AbstractChain{C,R}, e::C) where {C,R} = throw(MethodError(getindex, (typeof(ch), C)))

"""
    push!(ch:: AbstractChain{C,R}, e::Pair{C,R})

Add element with coefficient `e` to the chain `ch`. If the element is already in the chain then the coefficient
value of `e` will be added to the chain element coefficient.

Use `setindex!` to reset an element coefficient in the chain.
"""
push!(ch:: AbstractChain{C,R}, e::Pair{C,R}) where {C,R} = throw(MethodError(push!, (typeof(ch),Tuple{C,R})))

"""
    iterate(ch::AbstractChain)

Return a element/coefficient iterator of the chain `ch`
"""
iterate(ch::AbstractChain, state...) = throw(MethodError(simplify, (typeof(ch),Any)))

"""
    simplify(ch::AbstractChain)

Returns a simplified (without 0 elements) copy of the chain `ch`.
"""
simplify(ch::AbstractChain) = throw(MethodError(simplify, (typeof(ch),)))

# Auxiliary methods

function show(io::IO, ch::AbstractChain)
    print(io, "$(dim(ch)): ")
    if iszero(ch)
        print(io, "0")
    else
        for (i,(id,v)) in enumerate(ch.cells)
            i != 1 && print(io, " + ")
            print(io, "$v[0x$(lpad(string(id,base=16),16,'0'))]")
        end
    end
end

iszero(ch::AbstractChain) = length(ch) == 0
keytype(ch::AbstractChain{C,R}) where {C,R} = C
valtype(ch::AbstractChain{C,R}) where {C,R} = R
push!(ch::AbstractChain{C,R}, e::C, c::R) where {C,R} = push!(ch, e=>c)
push!(ch::AbstractChain{C,R}, e::Tuple{C,R}) where {C,R} = push!(ch, e[1]=>e[2])

"""
    map(f, ch::AbstractChain)

Apply function `f` to every coefficient of the chain `ch`.
"""
map(f, ch::AbstractChain) = map!(f, copy(ch), ch)

function map!(f, dest::AbstractChain, src::AbstractChain)
    for (k,v) in src
        dest[k] = f(v)
    end
    return dest
end

"""
    in(e, ch::AbstractChain)

Check if the element `e` in the chain `ch`.
"""
in(e::C, ch::AbstractChain{C,R}) where {C,R} = in(e, keys(ch))

"""
    append!(a::AbstractChain, b::AbstractChain)

Perform mutable addition of the chain `b` to the chain `a`.
"""
function append!(a::AbstractChain, b::AbstractChain)
    for (k,v) in b
        if k ∈ a
            a[k] += v
        else
            push!(a, k=>v)
        end
    end
    return a
end

"""
    +(ch::AbstractChain{C,R}, e::Tuple{C,R})

Perform non-mutable addition of the element with a coefficient `e` to the chain `ch`.
"""
(+)(ch::AbstractChain{C,R}, e::Tuple{C,R}) where {C,R} = push!(copy(ch), e)

"""
    +(a::AbstractChain, b::AbstractChain)

Perform non-mutable addition of the chain `b` to the chain `a`.
"""
(+)(a::AbstractChain, b::AbstractChain) = append!(copy(a), b)

"""
    +(a::AbstractChain, b::AbstractChain)

Perform non-mutable addition of the chain `b` to the chain `a`.
"""
function (-)(a::AbstractChain, b::AbstractChain)
    c = copy(a)
    for (k,v) in b
        push!(c, k => -v)
    end
    return c
end

"""
    +(ch::AbstractChain{C,R}, r::R)

Perform non-mutable addition of the coefficient `r` to every element of the chain `ch`.
"""
(+)(ch::AbstractChain{C,R}, r::R) where {C,R} = map(e->e+r, ch)
(+)(r::R, ch::AbstractChain{C,R}) where {C,R} = ch+r

"""
    *(ch::AbstractChain{C,R}, r::R)

Perform non-mutable addition of the coefficient `r` to every element of the chain `ch`.
"""
(*)(ch::AbstractChain{C,R}, r::R) where {C,R} = map(e->e*r, ch)
(*)(r::R, ch::AbstractChain{C,R}) where {C,R} = ch*r
