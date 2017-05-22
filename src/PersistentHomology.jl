module PersistentHomology

import Base: ==, -, +, *

export AbstractCell, Simplex, Cube,
       dim, faces, volume,

       Chain,
       setdim!, simplify

include("cells.jl")
include("chain.jl")

end # module
