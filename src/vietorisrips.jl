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
        # σ ← τ ∪ {v}
        first(v[:values]) in τvals && continue
        σ = Simplex(τvals..., v[:values]...)
        # σ ← τ ∪ {v}
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

"""Construction of Vietoris–Rips complex (an abstract simplicial complex)

This is a inefficient implementation (based on distance matrix) of the VRC algorithm from:

    Afra Zomorodian, Fast construction of the Vietoris-Rips complex, 2010
"""
function vietorisrips(X::AbstractMatrix, ɛ::Float64, weights=false;
                      method=:incremental)
    d, n = size(X)

    # add zero-dimensional simplices to complex
    splxs = map(Simplex, 1:n)
    cplx = SimplicialComplex(splxs...)

    # calculate cross distances for each point in dataset
    dist = Distances.pairwise(Distances.Euclidean(), X)

    # build 1-skeleton (neighborhood graph)
    E = spzeros(Bool, n, n) # adjacency matrix
    for σ in get(cells(cplx, 0))
        u = σ[:values]
        edges = find(d->d <= ɛ && d > 0., dist[:,u])
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
    kmax = maximum(mapslices(c->count(d->d <= ɛ, c), dist, 1))

    # perform VR expansion
    if method == :incremental
        incremental!(cplx, kmax, E)
    elseif method == :inductive
        inductive!(cplx, kmax, E)
    else
        throw(ArgumentError("Invalid method name $(method)"))
    end

    # calculate weights
    if !weights
        w = nothing
    else
        w = Dict{Int, Vector{Float64}}()
        w[0] = zeros(size(cplx,0))
        w[1] = zeros(size(cplx,1))
        for e in get(cells(cplx,1))
            w[1][e[:index]] = dist[e[:values]...]
        end
        for k in 2:dim(cplx)
            w[k] = zeros(size(cplx,k))
            for σ in get(cells(cplx,k))
                w[k][σ[:index]] = weight(σ, w, cplx)
            end
        end
    end

    return cplx, w
end
