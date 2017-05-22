module PersistentHomology

import Base: ==, -, +, *

export AbstractCell,
       dim, faces, volume,
       Simplex, Cube,

       AbstractChain,
       setdim!, simplify,
       Chain,

       AbstractComplex,
       boundary, coboundary, celltype, cells, boundary_matrix,
       SimplicialComplex,
       addsimplex, addsimplex!


include("cells.jl")
include("chain.jl")
include("complex.jl")
include("simplicial.jl")

end # module
