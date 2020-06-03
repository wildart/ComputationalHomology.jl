# Computational Homology

[![Build Status](https://travis-ci.org/wildart/ComputationalHomology.jl.svg?branch=master)](https://travis-ci.org/wildart/ComputationalHomology.jl)
[![Coverage Status](https://coveralls.io/repos/wildart/ComputationalHomology.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/wildart/ComputationalHomology.jl?branch=master)

This package provides various computational homology tools for cellular complexes.


## Installation
This package requires [SmithNormalForm.jl](https://github.com/wildart/SmithNormalForm.jl/) package for work.

```
pkg> add https://github.com/wildart/SmithNormalForm.jl.git#0.2.1
pkg> add https://github.com/wildart/ComputationalHomology.jl.git#0.2.0
```

For Julia 1.1+, add [BoffinStuff](https://github.com/wildart/BoffinStuff.git) registry in package manager, and proceed installation:

```
pkg> registry add https://github.com/wildart/BoffinStuff.git
pkg> add ComputationalHomology
```

## Features

- Cells

    - Simplex
    - Cube
    - CW

- Chains for specified PID

- Complexes

    - Simplicial
    - CW

- Filtrations

- Constructions

    - Čech
    - Vietoris-Rips
    - Witness

- Homology

    - Simplicial
    - Persistent

- Persistence

    - Barcodes / Diagrams
    - Persistence Landscape
    - Persistence Image
    - Distances
        - Wasserstein


## Example
```julia

julia> using ComputationalHomology

julia> X = rand(3,10); # generate dataset

julia> cplx, w = vietorisrips(X, 0.4, true) # generate Vietoris-Rips (VR) complex
(SimplicialComplex((10, 12, 4)), Dict(0=>[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],2=>[0.338893, 0.28014, 0.384243, 0.380966],1=>[0.338893, 0.321811, 0.304665, 0.310862, 0.27196, 0.28014, 0.366947, 0.380966, 0.191768, 0.384243, 0.359153, 0.365016]))

julia> flt = filtration(cplx, w) # construct filtration complex from VR complex
Filtration(SimplicialComplex((10, 12, 4)))

julia> ph = persistenthomology(flt, TwistReduction) # create persistent homology object with specific computation method
PersistentHomology[Filtration(SimplicialComplex((10, 12, 4))) with ComputationalHomology.TwistReduction]

julia> group(ph, 0) # calculate 0-homology group
2

julia> group(ph, 1) # calculate 1-homology group
3
```

## TODO
- [ ] Distances for persistance diagrams
- [ ] Landscape standard deviation
