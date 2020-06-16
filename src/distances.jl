rot(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

convert(::Type{Matrix}, dgm::PersistenceDiagram; skipinf=true) =
    hcat(( [first(i), last(i)] for i in dgm if !skipinf || !isinf(birth(i)) )...)

"""
    wasserstein(dgm1, dgm2)

Calculate Wasserstein distance between persistent diagrams `dgm1` and `dgm2`.
"""
function wasserstein(dgm1::PersistenceDiagram{T}, dgm2::PersistenceDiagram{T}) where {T}
    m1 = convert(Matrix, dgm1)
    m2 = convert(Matrix, dgm2)
    m, n = min(size(m1)...), min(size(m2)...)
    if m == 0
        m = 1
        m1 = zeros(T, 2, m)
    end
    if n == 0
        n = 1
        m2 = zeros(T, 2, n)
    end

    CSM = fill(0.0, n+m, n+m)
    ul = view(CSM,1:m,1:n)
    pairwise!(ul, Euclidean(), m1, m2, dims=2)
    maxd = maximum(ul)

    ur = @view CSM[1:m,(n+1):m+n]
    ur[:,:] .= maxd.*ones(m,m)
    view(ur, diagind(ur,0)) .= (m1'*rot(π/4))[:,2]

    ll = @view CSM[(m+1):m+n,1:n]
    ll[:,:] .= maxd.*ones(n,n)
    view(ll, diagind(ll,0)) .= (m2'*rot(π/4))[:,2]

    # Run the hungarian algorithm
    _, dist = hungarian(CSM)

    return dist
end

"""
    bottleneck(dgm1, dgm2)

Calculate bottleneck distance between persistent diagrams `dgm1` and `dgm2`.
"""
function bottleneck(dgm1::PersistenceDiagram{T}, dgm2::PersistenceDiagram{T}) where {T}
    return zero(T)
end
