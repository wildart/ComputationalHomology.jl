"""Bitmap Cubical Complex

A bitmap cubical complex is a structured complex representation of the collection of cubes.
The cubes store filtration values, thus this the cubical complexes is a filtered cubical complex.
"""
mutable struct BitmapComplex{T<:Real} <: AbstractComplex
    dims::Vector{Int}
    mult::Vector{Int}
    cells::AbstractVector{T} # cubes encoded in contiguous array
end

function BitmapComplex(::Type{T}, sizes::Vector{Int}=zeros(Int,0)) where {T<:Real}
    m = length(sizes) == 0 ? 0 : 1
    mlt = Int[]
    for s in sizes
        push!(mlt, m)
        m *= 2*s+1
    end
    BitmapComplex{T}(sizes, mlt, fill(typemax(T), m))
end
BitmapComplex() = BitmapComplex(Float64)

function BitmapComplex(sizes::Vector{Int}, cells::AbstractVector{T}) where {T<:Real}
    cplx = BitmapComplex(T, sizes)
    process = topcubes(cplx)
    for (i,v) in zip(process, cells)
        cplx.cells[i] = v
    end
    idx = BitVector(zeros(length(cplx.cells)))
    while length(process) > 0
        i = popfirst!(process)
        for j in faces(cplx, i)
            if cplx.cells[j] > cplx.cells[i]
                cplx.cells[j] = cplx.cells[i]
            end
            if !idx[j]
                push!(process, j)
                idx[j] = true
            end
        end
    end
    return cplx
end

BitmapComplex(img::AbstractMatrix{T}) where {T<:Real} = BitmapComplex([size(img)...], vec(rotr90(img)))

function topcubes(cplx::BitmapComplex{T}) where {T<:Real}
    idxs = Int[]
    pos=ones(Int, dim(cplx))
    while pos[1] <= cplx.dims[1]
        idx = cellindex(cplx, pos...)
        push!(idxs, idx)
        # println("pos: $pos, idx: $idx")
        pos[1]+=1
        i = findfirst(x->x>0, pos .- cplx.dims)
        if i !== nothing
            pos[i]=1
            pos[i+1]+=1
        end
        last(pos) > last(cplx.dims) && break
    end
    return idxs
end

cellindex(cplx::BitmapComplex, idxs...) = sum((2*i-1)*m for (i, m) in zip(idxs, cplx.mult))+1

function celldim(cplx::BitmapComplex, ci::Int)
    d = 0
    ci-=1
    for m in reverse(cplx.mult)
        pos = div(ci, m)
        if (pos % 2 == 1)
            d+=1
        end
        ci = mod(ci, m)
    end
    return d
end

function faces(cplx::BitmapComplex, ci::IX) where {IX<:Integer}
    idxs = IX[]
    d = 0
    c=ci-1
    for m in reverse(cplx.mult)
        pos = div(c, m)
        if (pos % 2 == 1)
            ms = d % 2 == 1 ? m : -m
            push!(idxs, ci + ms)
            push!(idxs, ci - ms)
            d+=1
        end
        c = mod(c, m)
    end
    return idxs
end

function cofaces(cplx::BitmapComplex, ci::IX) where {IX<:Integer}
    ret = IX[]
    c=ci-1
    p = foldl((p,m)-> (push!(p[1], div(p[2],m)+1), mod(p[2],m)),
              reverse(cplx.mult), init=(Int[],c)) |> first
    d = dim(cplx)+1
    for (i,m) in enumerate(reverse(cplx.mult))
        if div(c, m) % 2 == 0
            if ci > m+1 && p[i] != 1
                push!(ret, ci - m)
            end
            if (ci + m < length(cplx)) && (p[i] != 2 * cplx.dims[d-i]+1)
                push!(ret, ci + m)
            end
        end
        c = mod(c, m)
    end
    return ret
end

show(io::IO, cplx::BitmapComplex) = print(io, "BitmapComplex($(size(cplx)))")
eltype(cplx::BitmapComplex{T}) where {T} = T
dim(cplx::BitmapComplex) = length(cplx.dims)
length(cplx::BitmapComplex) = length(cplx.cells)
function size(cplx::BitmapComplex)
    ii = [celldim(cplx, i) for i in 1:length(cplx)]
    ntuple(d->count(i->i==d-1, ii), dim(cplx)+1)
end

# Filtration

function filtration(cplx::BitmapComplex{T}) where {T<:Real}
    U = promote_type(Float32, T)
    idx = Vector{Tuple{Int,UInt,U}}()
    i = 1
    for (ci,cv) in enumerate(cplx.cells)
        d = celldim(cplx, ci)
        push!(idx, (d, ci, cv))
    end
    sort!(idx, by=x->(x[3], x[1])) # sort by dimension & filtration value
    return Filtration(cplx, idx, U(Inf))
end

# Miscellaneous

function read(io::IO, ::Type{BitmapComplex}, ::Type{T}) where {T<:Real}
    d = parse(Int, readline(io))
    dims = Int[]
    for i in 1:d
        push!(dims, parse(Int, readline(io)))
    end
    n  = prod(dims)
    data = zeros(T, n)
    for i in 1:n
        data[i] =  parse(T, readline(io))
    end
    return BitmapComplex(dims, data)
end
