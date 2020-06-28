module ComputationalHomology

using LinearAlgebra
using Hungarian: hungarian
using SmithNormalForm: smith
using Distances: Metric, Euclidean, pairwise, colwise, pairwise!
using Random: shuffle!, randperm
using BoundingSphere: boundingsphere
using SparseArrays: SparseMatrixCSC, spzeros, findnz

import Base: ==, -, +, *, union, keys, values, hash, first, last, isless, show, length,
             eltype, valtype, getindex, setindex!, size, iterate, push!, append!, in, similar,
             read, write, vec, complex, iszero, convert, reduce, keytype, copy, map, map!,
             minimum, maximum, parse, isempty
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
       betti, euler,

       Filtration, filtration, order, simplices, similarity,

       AbstractPersistenceReduction,
       StandardReduction, TwistReduction,
       pairs, reduce!,
       PersistentHomology, persistenthomology,
       PersistentCocycleReduction, persistentcohomology, PersistentCocycleReduction,

       AbstractInterval, birth, death,
       Interval, AnnotatedInterval,
       PersistenceDiagram, diagram,
       Landscape, landscape, mean,
       PersistentImage,
       wasserstein

include("abstractchain.jl")
include("chains.jl")
include("abstractcell.jl")
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
include("iterators.jl")

@deprecate intervals(d, ps...) diagram(d, ps...)
@deprecate celltype(cplx) eltype(cplx)

end # module
