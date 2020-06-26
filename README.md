# Computational Homology

This package provides various computational homology tools for cellular complexes.

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url] [![][coverage-img]][coverage-url] |


## Installation

For Julia 1.1+, add [BoffinStuff](https://github.com/wildart/BoffinStuff.git) registry in package manager, and proceed with installation:

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

julia> cplx, w = vietorisrips(X, 0.4) # generate Vietoris-Rips (VR) complex
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


[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://wildart.github.io/ComputationalHomology.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://wildart.github.io/ComputationalHomology.jl/stable

[travis-img]: https://travis-ci.org/wildart/ComputationalHomology.jl.svg?branch=master
[travis-url]: https://travis-ci.org/wildart/ComputationalHomology.jl

[coverage-img]: https://img.shields.io/coveralls/wildart/ComputationalHomology.jl.svg
[coverage-url]: https://coveralls.io/r/wildart/ComputationalHomology.jl?branch=master

[issues-url]: https://github.com/wildart/ComputationalHomology.jl/issues
