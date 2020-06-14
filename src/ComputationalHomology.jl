module ComputationalHomology

using LinearAlgebra
using Hungarian: hungarian
using SmithNormalForm: smith
using Distances: Euclidean, pairwise, colwise, pairwise!
using Random: shuffle!, randperm
using BoundingSphere: boundingsphere
using SparseArrays: SparseMatrixCSC, spzeros, findnz

import Base: ==, -, +, *, union, values, hash, first, last, isless, show, length,
             eltype, valtype, getindex, size, iterate, push!, append!, in, read, write,
             vec, complex, iszero, convert, reduce
import Base.Iterators: pairs
import SparseArrays: sparse
import Statistics: mean
import LinearAlgebra: diag

export AbstractCell,
       dim, faces, volume,
       AbstractSimplex,
       values, vertices,
       Simplex, Cube, Cell,

       AbstractChain,
       Chain, simplify,

       AbstractComplex,
       boundary, coboundary, cells,
       SimplicialComplex, addsimplex!, addsimplices!,
       CWComplex,
       vietorisrips, witness, cech, čech,

       AbstractHomology, grouptype, group,
       Homology, homology, withgenerators, generators,

       Filtration, filtration, order, simplices, similarity,

       AbstractPersistenceReduction,
       StandardReduction, TwistReduction,
       pairs, reduce!,
       PersistentHomology, persistenthomology,

       Interval, birth, death,
       PersistenceDiagram, diagram,
       Landscape, landscape, mean,
       PersistentImage,
       wasserstein

include("abstractcell.jl")
include("chain.jl")
include("simplex.jl")
include("cube.jl")
include("cwcell.jl")
include("complex.jl")
include("simplicialcomplex.jl")
include("cwcomplex.jl")
include("constructions.jl")
include("filtration.jl")
include("homology.jl")
include("intervals.jl")
include("persistence.jl")
include("landscape.jl")
include("pimage.jl")
include("examples.jl")
include("distances.jl")

@deprecate intervals(d, ps...) diagram(d, ps...)
@deprecate celltype(cplx) eltype(cplx)

end # module
