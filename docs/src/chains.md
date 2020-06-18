# Chains

```@docs
AbstractChain
Chain{IX<:Integer, R}
ComputationalHomology.EmptyChain
```

Every chain implementation need to implement following interface functions.

```@docs
dim(::AbstractChain)
length(::AbstractChain)
keys(::AbstractChain)
values(::AbstractChain)
copy(::AbstractChain)
setindex!(::AbstractChain{C,R}, ::C, ::R) where {C,R}
getindex(::AbstractChain{C,R}, ::C) where {C,R}
push!(::AbstractChain{C,R}, ::Pair{C,R}) where {C,R}
iterate(::AbstractChain)
simplify(::AbstractChain)
```

There are additional function available after implementing the above interface.

```@docs
iszero(::AbstractChain)
keytype(::AbstractChain)
valtype(::AbstractChain)
push!(::AbstractChain{C,R}, ::C, ::R) where {C,R}
push!(::AbstractChain{C,R}, ::Tuple{C,R}) where {C,R}
map(f, ch::AbstractChain)
map!(f, ::AbstractChain, ::AbstractChain)
in(::C, ::AbstractChain{C,R}) where {C,R}
append!(::AbstractChain, ::AbstractChain)
(+)(::AbstractChain{C,R}, ::Tuple{C,R}) where {C,R}
(+)(::AbstractChain{C,R}, ::R) where {C,R}
```
