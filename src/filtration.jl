"""Filtration of an abstract complex

We call this sequence of complexes the **filtration** of `f` and
think of it as a construction by adding chunks of simplices at a time `t::FI`.
∅ = K0 ⊆ K1 ⊆ . . . ⊆ Kn = K.
"""
type Filtration{C<:AbstractComplex, FI}
    # underlying abstract cell complex
    complex::C
    # total order of simplexes defined by corresponding values of type FI
    index::Dict{FI,Vector{Tuple{Int,Int}}} # filtration_value => (cell dimension, cell index)
end
Base.show(io::IO, flt::Filtration) = print(io, "Filtration($(flt.complex))")
Base.length(flt::Filtration) = length(flt.index)

#
# Constructors
#

Filtration{C <: AbstractComplex, FI}(::Type{C}, ::Type{FI}) = Filtration(C(), Dict{FI,Vector{Tuple{Int,Int}}}())

function filtration{C <: AbstractComplex, FI}(cplx::C, ::Type{FI})
    idx = Dict{FI,Vector{Tuple{Int,Int}}}()
    i = one(FI)
    for d in 0:dim(cplx)
        for c in get(cells(cplx, d))
            idx[i] = [(dim(c), c[:index])]
            i += one(FI)
        end
    end
    return Filtration(cplx, idx)
end

"""Construct filtration from a cell complex and a complex weight function"""
function filtration{C <: AbstractComplex, FI}(cplx::C, w::Dict{Int,Vector{FI}})
    idx = Dict{FI,Vector{Tuple{Int,Int}}}()
    for d in 0:dim(cplx)
        for c in get(cells(cplx, d))
            ci = c[:index]
            fltval = w[d][ci]
            !haskey(idx, fltval) && setindex!(idx, FI[], fltval)
            push!(idx[fltval], (dim(c), c[:index]))
        end
    end
    return Filtration(cplx, idx)
end

function Base.push!{C<:AbstractComplex, FI}(flt::Filtration{C,FI}, cl::AbstractCell, v::FI; recursive=false)
    @assert isa(cl, celltype(flt.complex)) "Complex $(flt.complex) does not accept $(typeof(cl))"
    !haskey(flt.index, v) && setindex!(flt.index, Tuple{Int,Int}[], v)
    cls = push!(flt.complex, cl, recursive=recursive)
    for c in cls
        push!(flt.index[v], (dim(c), c[:index]))
    end
    return flt
end

"""Generate a combined boundary matrix from the filtration `flt` for the persistent homology calculations."""
function boundary_matrix(flt::Filtration; reduced=true)
    makereduced = convert(Int, reduced)
    # filtration total order map (simplex dimension => its order in dimension) => total order in filtration
    total = Dict{Pair{Int,Int}, Int}()
    # initialize boundary matrix
    bm = map(i->IntSet(), 1:sum(size(flt.complex))+makereduced)
    # fill boundary matrix
    col = 1 + makereduced
    for fltval in sort!(collect(keys(flt.index)))
        for (d, ci) in sort(flt.index[fltval], lt=(x,y)->(x[1] <= y[1])) # sort by dimension
            total[d=>ci] = col
            if d > 0
                splx = get(flt.complex[ci, d])
                for face in faces(splx)
                    fi = flt.complex[face, d-1]
                    push!(bm[col], total[d-1=>fi])
                end
            elseif reduced
                push!(bm[col], 1)
            end
            col += 1
        end
    end
    return bm
end

function Base.sparse(∂::Vector{IntSet})
    m = length(∂)
    ret = spzeros(Int, m, m)
    for i in 1:m
        bm = ∂[i]
        for (l, j) in enumerate(bm)
            ret[j,i] = j # (-1)^((l-1)%2) # coefs require exact order of faces in provides simplex
        end
    end
    collect(ret)
    return ret
end

"Write a combined boundary matrix to a text file"
function writeboundarymatrix(io::IO, bm::Vector, zeroindex = true)
    for smplx in bm
        if length(smplx) == 0
            write(io, "0")
        else
            write(io, "$(length(smplx)-1)")
        end
        for i in smplx
            write(io, " $(zeroindex ? i-1: i)")
        end
        write(io, 0x0A)
    end
end
