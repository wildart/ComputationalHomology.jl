#=
Each complex should have collection of cells per dimension:
- cells::Dict{Int,Vector{C}} or Vector{Vector{C}}
=#
abstract type AbstractComplex end

function boundary(cplx::AbstractComplex, ch::Chain{R}) where {R}
    d = dim(ch)
    cc = Chain(d-1, R)
    for (coef,elem) in ch
        append!(cc, coef * boundary(cplx, elem, d, R))
    end
    return cc
end

function coboundary(cplx::AbstractComplex, ch::Chain{R}) where {R}
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
boundary(cplx::AbstractComplex, i::Int, d::Int, ::Type{R}) where {R} = throw(MethodError(boundary, (typeof(cplx),Int,Int,R)))
boundary(cplx::AbstractComplex, i::Int, d::Int) = boundary(cplx, i, d, Int)

"""Return a complex coboundary given element and dimension"""
coboundary(cplx::AbstractComplex, i::Int, d::Int, ::Type{R}) where {R} = throw(MethodError(coboundary, (typeof(cplx),Int,Int,R)))
coboundary(cplx::AbstractComplex, i::Int, d::Int) = coboundary(cplx, i, d, Int)

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

"""
    push!(cplx::AbstractComplex, cell::AbstractCell; recursive=false) -> Vector{AbstractCell}

Insert a `cell` to a complex `cplx`, and returns an array of inserted cell(s). If `recursive=true` is passed then all faces of the `cell` are also added to the complex `cplx`.
"""
Base.push!(cplx::AbstractComplex, c::AbstractCell; recursive=false) = throw(MethodError(push!, (typeof(cplx),typeof(c))))

#
# Public Methods
#
"""Return a size of cell collections per dimension"""
Base.size(cplx::AbstractComplex) = (map(length, cells(cplx))...,)
Base.length(cplx::AbstractComplex) = sum(size(cplx))

"""Return a size of the cell collection for dimension (0-based)"""
function Base.size(cplx::AbstractComplex, d::Int)
    (d < 0 || d > dim(cplx)) && return 0
    sz = size(cplx)
    d >= length(sz) && return 0
    return sz[d+1]
end

"""Return an index of the k-dimensional cell in the complex.

**Note:** If returned index is larger than the number of cells of this dimension, then the cell is not in complex and the returned index will be assigned to the cell when it's added to the complex.
"""
function Base.getindex(cplx::AbstractComplex, c::C, d::Int) where {C}
    @assert celltype(cplx) == C "Incorrect cell type"
    dcells = cells(cplx, dim(c))
    dcells === nothing && return 1
    cidx = findfirst(isequal(c), dcells)
    cidx === nothing && return size(cplx, d)+1
    return dcells[cidx][:index]
end
Base.getindex(cplx::AbstractComplex, c::C) where {C} =  cplx[c, dim(c)]

"""Return a `d`-dimensional cell given its index `idx`."""
function Base.getindex(cplx::AbstractComplex, idx::Int, d::Int)
    cs = cells(cplx, d)
    cs === nothing && return nothing
    return cs[idx]
end

"""Generate a boundary matrix from the cell complex of the dimension `d`."""
function boundary_matrix(::Type{R}, cplx::AbstractComplex, d::Int) where {R}
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
boundary_matrix(cplx::AbstractComplex, d::Int) = boundary_matrix(Int, cplx, d)

#
# Complex simplex iterator
#
Base.length(splxs::Simplices{C}) where C<:AbstractComplex =
    splxs.dim < 0 ? sum(size(splxs.itr)) : size(splxs.itr, splxs.dim)

Base.eltype(iter::Simplices{C}) where C<:AbstractComplex = celltype(iter.itr)

# State (total id, dim, dim id)
function Base.iterate(iter::Simplices{C}, (tid, d, did)=(0, iter.dim, 1)) where C<:AbstractComplex
    tid >= length(iter) && return nothing
    if d < 0
        d = 0
    end
    if did > size(iter.itr, d)
        d += 1
        did = 1
    end
    return iter.itr[did, d], (tid+1, d, did+1)
end
