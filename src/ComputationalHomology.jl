module ComputationalHomology

using LinearAlgebra
import SparseArrays: SparseMatrixCSC, spzeros, sparse, findnz
import SmithNormalForm: smith
import Distances
import Random
import BoundingSphere
import Statistics: mean

import Base: ==, -, +, *, union, values, hash, first, last, isless, show, length,
             eltype, getindex, size, iterate, push!, append!, in, read, write,
             vec, complex, iszero
import Base.Iterators: pairs

export AbstractCell,
       dim, faces, volume,
       AbstractSimplex,
       values, vertices,
       Simplex, Cube, Cell,

       AbstractChain,
       Chain, simplify,

       AbstractComplex,
       boundary, coboundary, celltype, cells,
       SimplicialComplex, addsimplex!, addsimplices!,
       CWComplex,
       vietorisrips, witness, cech, čech,

       AbstractHomology, grouptype, group,
       Homology, homology, withgenerators, generators,

       Filtration, filtration, order, simplices,

       AbstractPersistenceReduction,
       StandardReduction, TwistReduction,
       pairs, intervals, reduce!, reduce,
       PersistentHomology, persistenthomology,

       Interval, PersistentDiagram, birth, death,
       Landscape, landscape, mean,
       PersistentImage

include("cells.jl")
include("chain.jl")
include("complex.jl")
include("simplicial.jl")
include("cw.jl")
include("cube.jl")
include("constructions.jl")
include("filtration.jl")
include("homology.jl")
include("persistence.jl")
include("landscape.jl")
include("pimage.jl")
include("examples.jl")

end # module
