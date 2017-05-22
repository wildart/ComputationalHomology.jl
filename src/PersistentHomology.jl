module PersistentHomology

import Base: ==, -, +, *

export AbstractCell, Simplex, Cube,
       dim, faces, volume,

       AbstractChain, Chain,
       setdim!, simplify,

       AbstractComplex,
       boundary, coboundary, celltype, cells, boundary_matrix

include("cells.jl")
include("chain.jl")
include("complex.jl")

end # module
