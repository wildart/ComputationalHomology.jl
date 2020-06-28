var documenterSearchIndex = {"docs":
[{"location":"cells/#Cells-1","page":"Cells","title":"Cells","text":"","category":"section"},{"location":"cells/#","page":"Cells","title":"Cells","text":"AbstractCell\ndim(::AbstractCell)\nfaces(::AbstractCell)\nhash(::AbstractCell)\nunion(::C, ::C) where {C<:AbstractCell}\nvertices(::AbstractCell)\nboundary(::Type{R}, ::AbstractCell) where {R}","category":"page"},{"location":"cells/#ComputationalHomology.AbstractCell","page":"Cells","title":"ComputationalHomology.AbstractCell","text":"Abstract cell type\n\n\n\n\n\n","category":"type"},{"location":"cells/#ComputationalHomology.dim-Tuple{AbstractCell}","page":"Cells","title":"ComputationalHomology.dim","text":"Dimension of the cell\n\n\n\n\n\n","category":"method"},{"location":"cells/#ComputationalHomology.faces-Tuple{AbstractCell}","page":"Cells","title":"ComputationalHomology.faces","text":"faces(cell)\n\nGet a collection of cell faces.\n\n\n\n\n\n","category":"method"},{"location":"cells/#Base.hash-Tuple{AbstractCell}","page":"Cells","title":"Base.hash","text":"hash(cell)\n\nGet a unique identifier of the cell\n\n\n\n\n\n","category":"method"},{"location":"cells/#Base.union-Union{Tuple{C}, Tuple{C,C}} where C<:AbstractCell","page":"Cells","title":"Base.union","text":"union(u, v)\n\nCreate a new cell from a combination of vertices of cells u and v.\n\n\n\n\n\n","category":"method"},{"location":"cells/#ComputationalHomology.vertices-Tuple{AbstractCell}","page":"Cells","title":"ComputationalHomology.vertices","text":"vertices(cell)\n\nGet a collection of cell vertices.\n\n\n\n\n\n","category":"method"},{"location":"cells/#ComputationalHomology.boundary-Union{Tuple{R}, Tuple{Type{R},AbstractCell}} where R","page":"Cells","title":"ComputationalHomology.boundary","text":"boundary(cell)\n\nCalculate boundary chain of the cell.\n\n\n\n\n\n","category":"method"},{"location":"cells/#Simplex-1","page":"Cells","title":"Simplex","text":"","category":"section"},{"location":"cells/#","page":"Cells","title":"Cells","text":"AbstractSimplex\nSimplex\nvalues(::AbstractSimplex)","category":"page"},{"location":"cells/#ComputationalHomology.AbstractSimplex","page":"Cells","title":"ComputationalHomology.AbstractSimplex","text":"Abstract simplex type\n\n\n\n\n\n","category":"type"},{"location":"cells/#Base.values-Tuple{AbstractSimplex}","page":"Cells","title":"Base.values","text":"values(cell)\n\nGet an array of cell values.\n\n\n\n\n\n","category":"method"},{"location":"cells/#CW-Cell-1","page":"Cells","title":"CW Cell","text":"","category":"section"},{"location":"cells/#","page":"Cells","title":"Cells","text":"Cell","category":"page"},{"location":"cells/#ComputationalHomology.Cell","page":"Cells","title":"ComputationalHomology.Cell","text":"CW cell type\n\n\n\n\n\n","category":"type"},{"location":"cells/#Cube-1","page":"Cells","title":"Cube","text":"","category":"section"},{"location":"cells/#","page":"Cells","title":"Cells","text":"Cube","category":"page"},{"location":"chains/#Chains-1","page":"Chains","title":"Chains","text":"","category":"section"},{"location":"chains/#","page":"Chains","title":"Chains","text":"AbstractChain\nChain{IX<:Integer, R}\nComputationalHomology.EmptyChain","category":"page"},{"location":"chains/#ComputationalHomology.AbstractChain","page":"Chains","title":"ComputationalHomology.AbstractChain","text":"AbstractChain{C, R}\n\nA chain is a formal linear combination of C-type elements in R-module.\n\nGiven a C-type element set E, we can construct a free R-module M that has E  M as a basis M, the module M is the module of the formal linear combinations of elements of E, or free module over E, R^(E), s.t. given a finite subset X_1  X_n of E, a formal linear combination of X_1  X_n is\n\na_1 X_1 +  + a_n X_n\n\nwhere the a_i ∈ R.\n\n\n\n\n\n","category":"type"},{"location":"chains/#ComputationalHomology.Chain","page":"Chains","title":"ComputationalHomology.Chain","text":"Chain{IX<:Integer, R}\n\nThe chain of integer elements for R-module.\n\nThis chain implementation uses a dictionary for storing elements and coefficients.\n\n\n\n\n\n","category":"type"},{"location":"chains/#ComputationalHomology.EmptyChain","page":"Chains","title":"ComputationalHomology.EmptyChain","text":"EmptyChain\n\nIt's a dummy chain class that does not store anything.\n\n\n\n\n\n","category":"type"},{"location":"chains/#","page":"Chains","title":"Chains","text":"Every chain implementation need to implement following interface functions.","category":"page"},{"location":"chains/#","page":"Chains","title":"Chains","text":"dim(::AbstractChain)\nlength(::AbstractChain)\nkeys(::AbstractChain)\nvalues(::AbstractChain)\ncopy(::AbstractChain)\nsetindex!(::AbstractChain{C,R}, ::C, ::R) where {C,R}\ngetindex(::AbstractChain{C,R}, ::C) where {C,R}\npush!(::AbstractChain{C,R}, ::Pair{C,R}) where {C,R}\niterate(::AbstractChain)\nsimplify(::AbstractChain)","category":"page"},{"location":"chains/#ComputationalHomology.dim-Tuple{AbstractChain}","page":"Chains","title":"ComputationalHomology.dim","text":"dim(ch::AbstractChain)\n\nReturns the dimension of the chain ch.\n\n\n\n\n\n","category":"method"},{"location":"chains/#Base.length-Tuple{AbstractChain}","page":"Chains","title":"Base.length","text":"length(ch::AbstractChain)\n\nReturns the length of the chain ch.\n\n\n\n\n\n","category":"method"},{"location":"chains/#Base.keys-Tuple{AbstractChain}","page":"Chains","title":"Base.keys","text":"keys(ch::AbstractChain)\n\nReturns elements of the chain ch.\n\n\n\n\n\n","category":"method"},{"location":"chains/#Base.values-Tuple{AbstractChain}","page":"Chains","title":"Base.values","text":"values(ch::AbstractChain)\n\nReturns coefficients of the chain ch.\n\n\n\n\n\n","category":"method"},{"location":"chains/#Base.copy-Tuple{AbstractChain}","page":"Chains","title":"Base.copy","text":"copy(ch::AbstractChain)\n\nReturns the copy of the chain ch.\n\n\n\n\n\n","category":"method"},{"location":"chains/#Base.getindex-Union{Tuple{R}, Tuple{C}, Tuple{AbstractChain{C,R},C}} where R where C","page":"Chains","title":"Base.getindex","text":"getindex(ch::AbstractChain{C,R}, e::C)\n\nReturns a coefficient of the element e of the chain ch.\n\n\n\n\n\n","category":"method"},{"location":"chains/#Base.push!-Union{Tuple{R}, Tuple{C}, Tuple{AbstractChain{C,R},Pair{C,R}}} where R where C","page":"Chains","title":"Base.push!","text":"push!(ch:: AbstractChain{C,R}, e::Pair{C,R})\n\nAdd element with coefficient e to the chain ch. If the element is already in the chain then the coefficient value of e will be added to the chain element coefficient.\n\nUse setindex! to reset an element coefficient in the chain.\n\n\n\n\n\n","category":"method"},{"location":"chains/#Base.iterate-Tuple{AbstractChain}","page":"Chains","title":"Base.iterate","text":"iterate(ch::AbstractChain)\n\nReturn a element/coefficient iterator of the chain ch\n\n\n\n\n\n","category":"method"},{"location":"chains/#ComputationalHomology.simplify-Tuple{AbstractChain}","page":"Chains","title":"ComputationalHomology.simplify","text":"simplify(ch::AbstractChain)\n\nReturns a simplified (without 0 elements) copy of the chain ch.\n\n\n\n\n\n","category":"method"},{"location":"chains/#","page":"Chains","title":"Chains","text":"There are additional function available after implementing the above interface.","category":"page"},{"location":"chains/#","page":"Chains","title":"Chains","text":"iszero(::AbstractChain)\nkeytype(::AbstractChain)\nvaltype(::AbstractChain)\npush!(::AbstractChain{C,R}, ::C, ::R) where {C,R}\npush!(::AbstractChain{C,R}, ::Tuple{C,R}) where {C,R}\nmap(f, ch::AbstractChain)\nmap!(f, ::AbstractChain, ::AbstractChain)\nin(::C, ::AbstractChain{C,R}) where {C,R}\nappend!(::AbstractChain, ::AbstractChain)\n(+)(::AbstractChain{C,R}, ::Tuple{C,R}) where {C,R}\n(+)(::AbstractChain{C,R}, ::R) where {C,R}","category":"page"},{"location":"chains/#Base.map-Tuple{Any,AbstractChain}","page":"Chains","title":"Base.map","text":"map(f, ch::AbstractChain)\n\nApply function f to every coefficient of the chain ch.\n\n\n\n\n\n","category":"method"},{"location":"chains/#Base.in-Union{Tuple{R}, Tuple{C}, Tuple{C,AbstractChain{C,R}}} where R where C","page":"Chains","title":"Base.in","text":"in(e, ch::AbstractChain)\n\nCheck if the element e in the chain ch.\n\n\n\n\n\n","category":"method"},{"location":"chains/#Base.append!-Tuple{AbstractChain,AbstractChain}","page":"Chains","title":"Base.append!","text":"append!(a::AbstractChain, b::AbstractChain)\n\nPerform mutable addition of the chain b to the chain a.\n\n\n\n\n\n","category":"method"},{"location":"chains/#Base.:+-Union{Tuple{R}, Tuple{C}, Tuple{AbstractChain{C,R},Tuple{C,R}}} where R where C","page":"Chains","title":"Base.:+","text":"+(ch::AbstractChain{C,R}, e::Tuple{C,R})\n\nPerform non-mutable addition of the element e with a coefficient to the chain ch.\n\n\n\n\n\n","category":"method"},{"location":"chains/#Base.:+-Union{Tuple{R}, Tuple{C}, Tuple{AbstractChain{C,R},R}} where R where C","page":"Chains","title":"Base.:+","text":"+(ch::AbstractChain{C,R}, r::R)\n\nPerform non-mutable addition of the coefficient r to every element of the chain ch.\n\n\n\n\n\n","category":"method"},{"location":"intervals/#Intervals-1","page":"Intervals","title":"Intervals","text":"","category":"section"},{"location":"intervals/#","page":"Intervals","title":"Intervals","text":"Any interval type need to implement following interface:","category":"page"},{"location":"intervals/#","page":"Intervals","title":"Intervals","text":"AbstractInterval\nbirth\ndeath","category":"page"},{"location":"intervals/#ComputationalHomology.AbstractInterval","page":"Intervals","title":"ComputationalHomology.AbstractInterval","text":"Abstract interval type\n\n\n\n\n\n","category":"type"},{"location":"intervals/#ComputationalHomology.birth","page":"Intervals","title":"ComputationalHomology.birth","text":"birth(i::AbstractInterval)\n\nReturn a birth value of the interval i.\n\n\n\n\n\n","category":"function"},{"location":"intervals/#ComputationalHomology.death","page":"Intervals","title":"ComputationalHomology.death","text":"death(i::AbstractInterval)\n\nReturn a death value of the interval i.\n\n\n\n\n\n","category":"function"},{"location":"intervals/#","page":"Intervals","title":"Intervals","text":"The persistence diagram is defined as","category":"page"},{"location":"intervals/#","page":"Intervals","title":"Intervals","text":"PersistenceDiagram{T} = AbstractVector{<:AbstractInterval{T}}","category":"page"},{"location":"intervals/#","page":"Intervals","title":"Intervals","text":"Following auxiliary functions are available for any interval instance derived from AbstractInterval:","category":"page"},{"location":"intervals/#","page":"Intervals","title":"Intervals","text":"birthx(::AbstractInterval)\ndeathx(::AbstractInterval)\npair(::AbstractInterval)\nisempty(::AbstractInterval)\nin(::T, ::AbstractInterval{T}) where {T<:AbstractFloat}\nisless(::AbstractInterval, ::AbstractInterval)\n(==)(::AbstractInterval, ::AbstractInterval)","category":"page"},{"location":"intervals/#","page":"Intervals","title":"Intervals","text":"Implemented interval types:","category":"page"},{"location":"intervals/#","page":"Intervals","title":"Intervals","text":"Interval\nAnnotatedInterval","category":"page"},{"location":"intervals/#ComputationalHomology.Interval","page":"Intervals","title":"ComputationalHomology.Interval","text":"Interval{T<:AbstractFloat} <: AbstractInterval{T}\n\nSimple implementation of the AbstractInterval type.\n\n\n\n\n\n","category":"type"},{"location":"intervals/#ComputationalHomology.AnnotatedInterval","page":"Intervals","title":"ComputationalHomology.AnnotatedInterval","text":"AnnotatedInterval{T<:AbstractFloat, C<:AbstractChain} <: AbstractInterval{T}\n\nInterval type annotated with a generator chain.\n\n\n\n\n\n","category":"type"},{"location":"#ComputationalHomology.jl-1","page":"Home","title":"ComputationalHomology.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"The package ComputationalHomology provides various computational homology tools for cellular complexes.","category":"page"},{"location":"#Getting-started-1","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"For Julia 1.1+, add BoffinStuff registry in package manager, and proceed with installation:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"pkg> registry add https://github.com/wildart/BoffinStuff.git\npkg> add ComputationalHomology","category":"page"},{"location":"#","page":"Home","title":"Home","text":"A simple example of computing the persistenthomology of the Vietoris–Rips complex.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"using ComputationalHomology\nX = rand(3,10); # generate dataset\nflt = filtration(vietorisrips(X, 0.4)...)\nph = persistenthomology(flt)\ngroup(ph, 0)","category":"page"}]
}
