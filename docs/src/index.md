# ComputationalHomology.jl

The package __ComputationalHomology__ provides various computational homology tools for cellular complexes.

## Getting started

For Julia 1.1+, add [BoffinStuff](https://github.com/wildart/BoffinStuff.git) registry in package manager, and proceed with installation:

```julia
pkg> registry add https://github.com/wildart/BoffinStuff.git
pkg> add ComputationalHomology
```

A simple example of computing the [`persistenthomology`](@ref) of the Vietoris–Rips complex.

```@repl
using ComputationalHomology
X = rand(3,10); # generate dataset
flt = filtration(vietorisrips(X, 0.4)...)
ph = persistenthomology(flt)
group(ph, 0)
```
