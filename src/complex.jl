#=
Each complex should have collection of cells per dimension:
- cells::Dict{Int,Vector{C}} or Vector{Vector{C}}
=#
abstract type AbstractComplex end

#
# AbstractComplex Public Interface
#
"""Return a complex boundary given element and dimension"""
boundary(cplx::AbstractComplex, i::Integer, d::Int, ::Type{PID}) where {PID} =
    throw(MethodError(boundary, (typeof(cplx),Int,Int,PID)))
boundary(cplx::AbstractComplex, i::Integer, d::Int) = boundary(cplx, i, d, Int)
boundary(cplx::AbstractComplex, splx::AbstractCell, ::Type{PID}) where {PID} =
    boundary(cplx, hash(splx), dim(splx), PID)
boundary(cplx::AbstractComplex, splx::AbstractCell) = boundary(cplx, hash(splx), dim(splx), Int)

"""Return a complex coboundary given element and dimension"""
coboundary(cplx::AbstractComplex, i::Integer, d::Int, ::Type{PID}) where {PID} =
    throw(MethodError(coboundary, (typeof(cplx),Int,Int,PID)))
coboundary(cplx::AbstractComplex, i::Integer, d::Int) = coboundary(cplx, i, d, Int)
coboundary(cplx::AbstractComplex, splx::AbstractCell, ::Type{PID}) where {PID} =
    coboundary(cplx, hash(splx), dim(splx), PID)
coboundary(cplx::AbstractComplex, splx::AbstractCell) = coboundary(cplx, hash(splx), dim(splx), Int)

"""Return a complex cell type"""
eltype(cplx::AbstractComplex) = throw(MethodError(eltype, (typeof(cplx),)))

"""Return a cell collection per dimension (increasing)"""
cells(cplx::AbstractComplex) = throw(MethodError(cells, (typeof(cplx),)))

"""Return a nullable cell collection per dimension (increasing)"""
cells(cplx::AbstractComplex, d::Int) = throw(MethodError(cells, (typeof(cplx),Int)))

"""
    push!(cplx::AbstractComplex, cell::AbstractCell; recursive=false) -> Vector{AbstractCell}

Insert a `cell` to a complex `cplx`, and returns an array of inserted cell(s). If `recursive=true` is passed then all faces of the `cell` are also added to the complex `cplx`.
"""
push!(cplx::AbstractComplex, c::AbstractCell; recursive=false) =
    throw(MethodError(push!, (typeof(cplx),typeof(c))))

#
# Public Methods
#
"""Return a number of cells per dimension"""
size(cplx::AbstractComplex) = (map(length, cells(cplx))...,)

"""Return a dimension of the complex"""
dim(cplx::AbstractComplex) = length(size(cplx))-1

"""Return a total number of cells in the complex"""
length(cplx::AbstractComplex) = sum(size(cplx))

"""Return a number of the cell in the complex of a dimension `d` (0-based)"""
function size(cplx::AbstractComplex, d::Int)
    sz = size(cplx)
    szlen = length(sz)
    (d < 0 || d >= szlen) && return 0
    return sz[d+1]
end

"""
    position(complex, index, dimenion)

Return a position of the cell in an order of cells of the same dimenion of the `complex` given its `index` and `dimenion`.
"""
function position(cplx::AbstractComplex, idx::Integer, d::Int)
    dcells = cells(cplx, d)
    length(dcells) == 0 && return 0
    cidx = findfirst(c->hash(c) == idx, dcells)
    return cidx === nothing ? 0 : cidx
end

"""
    position(complex, cell)

Return a position of the `cell` in an order of cells of the same dimenion of the `complex`.
"""
position(cplx::AbstractComplex, c::AbstractCell) = position(cplx, hash(c), dim(c))

"""
    cplx[idx, d]

Return a `d`-dimensional cell given its index `idx` and dimenion `d`.
"""
function getindex(cplx::AbstractComplex, idx::Integer, d::Int)
    cidx = position(cplx, idx, d)
    cidx == 0 && return nothing
    return cells(cplx, d)[cidx]
end

"""
    boundary(cplx, ch)

Return the chain `ch` boundary in the complex `cplx`.
"""
function boundary(cplx::AbstractComplex, ch::Chain{IX,R}) where {R, IX<:Integer}
    d = dim(ch)
    cc = Chain(d-1, IX, R)
    for (elem, coef) in ch
        append!(cc, coef * boundary(R, cplx[elem, d]))
    end
    return simplify(cc)
end

"""
    coboundary(cplx, ch)

Return the chain `ch` coboundary in the complex `cplx`.
"""
function coboundary(cplx::AbstractComplex, ch::Chain{IX,R}) where {R, IX<:Integer}
    d = dim(ch)
    cc = Chain(d+1, IX, R)
    # δₙ₋₁(c)(σ) = c(∂ₙ(σ))
    for (elem, coef) in ch
        append!(cc, coef * coboundary(cplx, elem, d, R))
    end
    return simplify(cc)
end

"""
    boundary(complex, d, PID)

Generate a boundary matrix from the cell `complex` of the dimension `d` in using coefficients of `PID`.
"""
function boundary(cplx::AbstractComplex, d::Int, ::Type{PID}) where {PID}
    csize = size(cplx)
    rows = d > 0 ? csize[d] : 0
    cols = d <= dim(cplx) ? csize[d+1] : 0
    bm = spzeros(PID, rows, cols)
    if d>=0 && d <= dim(cplx)
        for c in cells(cplx, d)
            i = position(cplx, hash(c), d)
            for (elem, coef) in boundary(PID, c)
                j = position(cplx, elem, d-1)
                bm[j, i] = coef
            end
        end
    end
    return bm
end

"""
    cell in cplx

Checks if the `cell` is in the complex `cplx`
"""
in(c::AbstractCell, cplx::AbstractComplex) = position(cplx, c) > 0

"""
    cochain(complex, d, coefficients)

Return a cochain of the dimension `d` for the `complex` with PID `coefficients`.
"""
function cochain(cplx::AbstractComplex, d::Int, coefs::Vector{PID}) where {PID}
    cs = cells(cplx, d)
    @assert length(cs) == length(coefs) "Number of coefficients must match number of cells of dimension $d"
    return Chain(d, coefs, map(c->hash(c), cs))
end
cochain(cplx::AbstractComplex, d::Int, coef::PID) where {PID} =
    cochain(cplx, d, fill(coef, size(cplx, d)))

"""
    adjacency_matrix(cplx, [T=Int])

Construct an adjacency matrix of type `T` from a 1-skeleton (1D subcomplex) of the complex `cplx`.
"""
function adjacency_matrix(cplx::AbstractComplex, ::Type{T}) where {T<:Integer}
    C0 = map(hash, cells(cplx, 0))
    N = length(C0)
    adj = spzeros(T,N,N)
    for c in cells(cplx, 1)
        i, j = map(h->findfirst(isequal(h), C0), vertices(c))
        adj[i, j] = adj[j, i] = one(T)
    end
    return adj
end
adjacency_matrix(cplx::AbstractComplex) = adjacency_matrix(cplx, Int)

"""
    showchain(cplx::AbstractComplex, ch::AbstractChain)

Display the chain `ch` with elements from the complex `cplx`.
"""
function showchain(cplx::AbstractComplex, ch::AbstractChain)
    d = dim(ch)
    R = valtype(ch)
    pos = one(R)
    print("[$d]: ")
    if iszero(ch)
        print("0")
    else
        for (i,(id,v)) in enumerate(ch.cells)
            val = abs(v)
            if sign(v) == pos
                i != 1 && print(" + ")
            else
                print(" - ")
            end
            splx = strip(repr("text/plain", cplx[id, d]), ['\"'])
            print("$val⋅$splx")
        end
    end
end
