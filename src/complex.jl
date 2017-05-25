#=
Each complex should have collection of cells per dimension:
- cells::Dict{Int,Vector{C}} or Vector{Vector{C}}
=#
abstract AbstractComplex

function boundary{R}(cplx::AbstractComplex, ch::Chain{R})
    d = dim(ch)
    cc = Chain(d-1,R)
    for (coef,elem) in ch
        append!(cc, coef * boundary(cplx, elem, d, R))
    end
    return cc
end

function coboundary{R}(cplx::AbstractComplex, ch::Chain{R})
    d = dim(ch)
    cc = Chain(d+1,R)
    for (coef,elem) in ch
        append!(cc, coef * coboundary(cplx, elem, d, R))
    end
    return cc
end

#
# AbstractComplex Public Interface
#
"""Return a complex boundary given element and dimension"""
boundary{M}(cplx::AbstractComplex, i::M, d::Int) = throw(MethodError(boundary, (typeof(cplx),M,Int)))

"""Return a complex coboundary given element and dimension"""
coboundary{M}(cplx::AbstractComplex, i::M, d::Int) = throw(MethodError(coboundary, (typeof(cplx),M,Int)))

"""Return a complex cell type"""
celltype(cplx::AbstractComplex) = throw(MethodError(celltype, (typeof(cplx),)))

"""Return a dimension of the complex"""
dim(cplx::AbstractComplex) = throw(MethodError(dim,(typeof(cplx),)))

"""Set dimension of the complex"""
setdim!(cplx::AbstractComplex, d::Int) = throw(MethodError(setdim!, (typeof(cplx),Int)))

"""Return a cell collection per dimension (increasing)"""
cells(cplx::AbstractComplex) = throw(MethodError(cells, (typeof(cplx),)))

"""Return a nullable cell collection per dimension (increasing)"""
cells(cplx::AbstractComplex, d::Int) = throw(MethodError(cells, (typeof(cplx),Int)))

#
# Public Methods
#
"""Return a size of cell collections per dimension"""
Base.size(cplx::AbstractComplex) = (map(length, cells(cplx))...)

"""Return a size of the cell collection for dimension (0-based)"""
function Base.size{C<:AbstractComplex}(cplx::C, d::Int)
    (d < 0 || d > dim(cplx)) && return 0
    sz = size(cplx)
    d >= length(sz) && return 0
    return sz[d+1]
end

"""Return an index of the k-dimensional cell in the complex.

**Note:** If returned index is larger than the number of cells of this dimension, then the cell is not in complex and the returned index will be assigned to the cell when it's added to the complex.
"""
function Base.getindex{C}(cplx::AbstractComplex, c::C, d::Int)
    @assert celltype(cplx) == C "Incorrect cell type"
    @assert dim(c) == d "Incorrect cell dimension"
    ndcells = cells(cplx, d)
    isnull(ndcells) && return 1
    dcells = get(ndcells)
    cidx = findfirst(dcells, c)
    cidx == 0 && return size(cplx, d)+1
    return dcells[cidx][:index]
end
Base.getindex{C}(cplx::AbstractComplex, c::C) =  cplx[c, dim(c)]

"""Return a `d`-dimensional cell given its index `idx`."""
function Base.getindex(cplx::AbstractComplex, idx::Int, d::Int)
    cs = cells(cplx, d)
    isnull(cs) ? Nullable{eltype(eltype(cs))}() : Nullable(get(cs)[idx])
end

"""Generate a boundary matrix from the cell complex of the dimension `d`."""
function boundary_matrix{R}(::Type{R}, cplx::AbstractComplex, d::Int)
    csize = size(cplx)
    rows = d > 0 ? csize[d] : 0
    cols = d <= dim(cplx) ? csize[d+1] : 0
    bm = spzeros(R, rows, cols)
    if d>=0 && d <= dim(cplx)
        for i in 1:csize[d+1]
            for (coef,elem) in boundary(cplx, i, d, R)
                bm[elem, i] = coef
            end
        end
    end
    return bm
end

# function Base.setindex!{C}(cplx::AbstractComplex, c::C, d::Int)
#     @assert celltype(cplx) == C "Incorrect cell type"
#     cidx = indexes(cplx)
#     length(cidx) <= d && resize!(cidx, d+1)
#     ccel = cells(cplx)
#     length(ccel) <= d && resize!(ccel, d+1)
#     length(size(cplx)) <= d && size!(cplx, 0, d+1)

#     if !haskey(cidx[d], c)
#         cidx[d][c] = length(cidx[d]) # index stet to last element
#         push!(cidx[d], c)
#         size!(cplx, length(cidx[d]), d)
#     end

#     return cplx
# end
