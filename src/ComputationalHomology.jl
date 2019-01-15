module ComputationalHomology

using LinearAlgebra
using SparseArrays
import SmithNormalForm: smith
import Distances
import Random
import BoundingSphere
import Statistics: mean

import Base: ==, -, +, *
import Base.Iterators: pairs

export AbstractCell,
       dim, faces, volume,
       Simplex, Cube,

       AbstractChain,
       setdim!, simplify,
       Chain,

       AbstractComplex,
       boundary, coboundary, celltype, cells, boundary_matrix,
       SimplicialComplex,
       vietorisrips, witness, cech, čech,

       AbstractHomology, grouptype, group,
       Homology, homology, withgenerators, generators,

       Filtration, filtration, order, simplices,

       AbstractPersistenceReduction,
       StandardReduction, TwistReduction,
       pairs, intervals, reduce!, reduce,
       PersistentHomology, persistenthomology,

       Interval, Landscape, landscape, mean

include("cells.jl")
include("chain.jl")
include("complex.jl")
include("simplicial.jl")
include("constructions.jl")
include("filtration.jl")
include("homology.jl")
include("persistence.jl")
include("landscape.jl")

end # module
