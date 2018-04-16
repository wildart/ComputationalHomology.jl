"""Find all neighbors of vertex `u`` within `G` that precede it in the given ordering"""
function lowernbrs(cplx, u, E)
    # get all edges of u
    uval = first(u[:values])
    uidx = u[:index]
    # select {v} s.t. u > v for all edges {u,v}
    !any(E[:,uval]) && return celltype(cplx)[]
    return filter(v->v[:index] < uidx && E[first(v[:values]),uval], get(cells(cplx,0)))
end

"""Incremental construction of VR complex from neighborhood graph"""
function inductive!(cplx, k, E)
    for i in 1:k
        cls = cells(cplx, i)
        isnull(cls) && continue
        for τ in get(cells(cplx, i))
            N = celltype(cplx)[]
            τvals = τ[:values]
            for (i, uval) in enumerate(τvals)
                uidx = cplx[Simplex(uval), 0]
                u = get(cplx[uidx, 0])
                N = i == 1 ? lowernbrs(cplx, u, E) : intersect(N, lowernbrs(cplx, u, E))
            end
            for v in N
                first(v[:values]) in τvals && continue
                σ = Simplex(τvals..., v[:values]...)
                cplx[σ, i] > size(cplx, i) && push!(cplx, σ)
            end
        end
    end
end

"""Add all simplices of the `k`-skeleton whose maximal vertex is `τ`"""
function addcofaces!(cplx, k, τ, N, E)
    # V ← V ∪ {τ}
    τdim = dim(τ)
    cplx[τ, τdim] > size(cplx, τdim) && push!(cplx, τ)
    # stop recursion
    dim(τ) >= k && return
    τvals = τ[:values]
    for v in N
        first(v[:values]) in τvals && continue
        # σ ← τ ∪ {v}
        σ = Simplex(τvals..., v[:values]...)
        M = intersect(N, lowernbrs(cplx, v, E))
        addcofaces!(cplx, k, σ, M, E)
    end
end

"""Incremental construction of VR complex from neighborhood graph"""
function incremental!(cplx, k, E)
    # construct VR complex
    for u in get(cells(cplx, 0))
        N = lowernbrs(cplx, u, E)
        addcofaces!(cplx, k, u, N, E)
    end
end

"""Weight function for VR complex"""
function weight(σ, w, cplx)
    k = dim(σ)
    k == 0 && return 0.
    k == 1 && return w[1][cplx[σ, 1]]
    return maximum([w[k-1][cplx[τ, k-1]] for τ in faces(σ)])
end

function expand(method, cplx, w, kmax, E)
    # perform nerve expansion
    if method == :incremental
        incremental!(cplx, kmax, E)
    elseif method == :inductive
        inductive!(cplx, kmax, E)
    else
        throw(ArgumentError("Invalid method name $(method)"))
    end

    if w !== nothing
        for k in 2:dim(cplx)
            w[k] = zeros(size(cplx,k))
            for σ in get(cells(cplx,k))
                w[k][σ[:index]] = weight(σ, w, cplx)
            end
        end
    end
    return
end

"""Construction of Vietoris–Rips complex (an abstract simplicial complex)

This is an inefficient implementation (based on distance matrix) of the VRC algorithm from:

    Afra Zomorodian, Fast construction of the Vietoris-Rips complex, 2010
"""
function vietorisrips(X::AbstractMatrix, ɛ::Float64, weights = true;
                      expansion=:incremental, distance = Distances.Euclidean(), maxdim = size(X,2)-1)
    d, n = size(X)

    # add zero-dimensional simplices to complex
    splxs = map(Simplex, 1:n)
    cplx = SimplicialComplex(splxs...)

    # calculate cross distances for each point in dataset
    D = Distances.pairwise(distance, X)

    # build 1-skeleton (neighborhood graph)
    E = spzeros(Bool, n, n) # adjacency matrix
    for σ in get(cells(cplx, 0))
        u = σ[:values]
        edges = find(d->d <= ɛ && d > 0., D[:,u])
        length(edges) == 0 && continue
        for i in edges
            v = get(cplx[i, 0])[:values]
            s = Simplex(u...,v...)
            if cplx[s] > size(cplx, 1)
                # add simplex to complex
                push!(cplx, s)
                # fill adjacency matrix
                E[s[:values], s[:values]] = true
            end
        end
    end

    # determine maximal dimension
    kmax = min(maxdim, maximum(mapslices(c->count(d->0.0<d≤ɛ, c), D, 1)))

    # calculate weights of nerve
    w = nothing
    if weights
        w = Dict{Int, Vector{Float64}}()
        w[0] = zeros(size(cplx,0))
        w[1] = zeros(size(cplx,1))
        for e in get(cells(cplx,1))
            w[1][e[:index]] = D[e[:values]...]
        end
    end

    # perform expansion
    expand(expansion, cplx, w, kmax, E)

    return cplx, w
end

"""
Landmark selection for witness complex
"""
function landmarks(X::AbstractMatrix, l::Int; method = :minmax,
                   distance = Distances.Euclidean(), firstpoint = 0)
    d, n = size(X)
    if method == :random
        return randperm(n)[1:l]
    elseif method == :minmax
        idxs = firstpoint > 0 ? [firstpoint] : rand(1:n, 1)
        Z = view(X, :, idxs)
        for i in 2:l
            maxv = -Inf
            maxvidx = 0
            for j in 1:n
                if j ∉ idxs
                    v, _ = findmin(Distances.colwise(distance, X[:,j], Z))
					if v > maxv
					    maxv = v
                        maxvidx = j
                    end
                end
            end
            push!(idxs, maxvidx)
        end
        return idxs
    else
        throw(ArgumentError("Invalid method name $(method)"))
    end
end

"""Construction of witness complex (an abstract simplicial complex)

This is an implementation of the witness construction algorithm from:

    Vin de Silva and Gunnar Carlsson, Topological estimation using witness complexes, 2004

For parameter ν = 0, 1, 2 determines size of landmark radius.
- If ν = 0, then for i = 1, 2, . . . , N define m_i = 0.
- If ν > 0, then for i = 1, 2, . . . , N define m_i to be the ν-th smallest entry of the i-th column of D.
"""
function witness(X::AbstractMatrix, l::Int, R::Float64, weights = true;
                 landmark = :minmax, distance = Distances.Euclidean(),
                 expansion = :incremental, ν::Int = 2, maxdim = size(X,2)-1,
                 firstpoint = 0)

    # get landmarks
    L = landmarks(X, l, method = landmark, distance = distance, firstpoint=firstpoint)

    # get distances to landmarks
    D = Distances.pairwise(distance, X[:,L], X)

    # simplexes contain indexes to landamarks
    cplx, w = witness(D, R, weights, expansion=expansion, ν=ν, maxdim=maxdim, firstpoint=firstpoint)

    return cplx, w, L
end

"""Genenerate witness complex from the distance matrix with simplexes constructed from landmark indexes"""
function witness(D::AbstractMatrix, R::Float64, weights = true;
                 expansion = :incremental, ν::Int = 2, maxdim = size(D,1)-1,
                 firstpoint = 0)

    @assert ν < 3 "ν only can take values of 0, 1 or 2"
    @assert maxdim > 0 "Maximum dimension should be more then 0"

    l, n = size(D)

    # get landmark radius extension m_i
    m = zeros(eltype(D), n)
    if ν > 0
        idxs = zeros(Int, l)
        for i in 1:n
            sortperm!(idxs, D[:,i])
            m[i] = D[idxs[ν], i]
        end
    end

    # construct 1-seleton edges
    splxs = map(Simplex, 1:l)
    cplx = SimplicialComplex(splxs...)
    w = nothing
    if weights
        w = Dict{Int, Vector{Float64}}()
        w[0] = zeros(size(cplx,0))
        w[1] = zeros(0)
    end

    # the edge σ = [ab] belongs to W(D; R, ν) iff there exists a witness i ∈ {1, 2, ..., N} such that:
    # max(D(a, i), D(b, i)) ≤ R + m_i
    E = spzeros(Bool, l, l) # adjacency matrix
    for i in 1:l
        for j in i+1:l
            e = Inf
            @inbounds for k in 1:n
                d = max(D[i, k], D[j, k])
                d = d < m[k] ? 0.0 : d - m[k]
                if d < e
                    e = d
                end
            end

            if e ≤ R
                push!(cplx, Simplex(i, j))
                weights && push!(w[1], e)
                E[i, j] = true
            end
        end
    end

    # determine maximal dimension
    kmax = min(maxdim, maximum(mapslices(c->count(d->0.0<d≤R, c), D, 1)))

    # perform expansion
    expand(expansion, cplx, w, kmax, E)

    return cplx, w
end
