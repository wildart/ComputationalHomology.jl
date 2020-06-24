# Cells

```@docs
AbstractCell
dim(::AbstractCell)
faces(::AbstractCell)
hash(::AbstractCell)
union(::C, ::C) where {C<:AbstractCell}
vertices(::AbstractCell)
boundary(::Type{R}, ::AbstractCell) where {R}
```

## Simplex

```@docs
AbstractSimplex
Simplex
values(::AbstractSimplex)
```

## CW Cell

```@docs
Cell
```

## Cube

```@docs
Cube
```
