# Intervals

Any interval type need to implement following interface:

```@docs
AbstractInterval
birth
death
```

The persistence diagram is defined as

```julia
PersistenceDiagram{T} = AbstractVector{<:AbstractInterval{T}}
```


Following auxiliary functions are available for any interval instance derived from
`AbstractInterval`:

```@docs
birthx(::AbstractInterval)
deathx(::AbstractInterval)
pair(::AbstractInterval)
isempty(::AbstractInterval)
in(::T, ::AbstractInterval{T}) where {T<:AbstractFloat}
isless(::AbstractInterval, ::AbstractInterval)
(==)(::AbstractInterval, ::AbstractInterval)
```

Implemented interval types:

```@docs
Interval
AnnotatedInterval
```
