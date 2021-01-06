const Orbits = [(1,0), (2,0), (0,1), (1,1), (0,2)]

"""
    flatten(pdm::Dict{Int, <:PersistenceDiagram}) -> PersistenceDiagram

Flatten a multi-dimensional pd into a barcode.
"""
flatten(pdm::Dict{Int, PD}) where {T<:Real, PD<:PersistenceDiagram{T}} =
    (collect∘Iterators.flatten∘values)(pdm)

"""
    maxsize(pds) -> Int

Return the number of bars in the largest barcode from the collection of barcodes `pds`
"""
maxsize(pds::AbstractVector{PD}) where {T<:Real, PD<:PersistenceDiagram{T}} =
    maximum(map(length, pds))

tropicdim(n) = n+(div(n*(n+1),2)) # codomain dimension

function tropicstats(pds::AbstractVector{PD}) where {T<:Real, PD<:PersistenceDiagram{T}}
    n = ComputationalHomology.maxsize(pds)
    d = ComputationalHomology.tropicdim(n)
    maxval = filter!(!isinf, map(i->max(birth(i),death(i)), Iterators.flatten(pds))) |> maximum
    m =  map(i->birth(i)/-ComputationalHomology.birthx(i), Iterators.flatten(pds)) |> maximum
    (d, n, ceil(Int, m), maxval)
end

"""
Caluclate 2-symmetric tropical rational polynomial for the persistence diagram
"""
function tropicp(pd::PersistenceDiagram{T},m,l,p; maxval=T(Inf)) where {T<:Real}
    B = Tuple{T,T}[]
    for i in pd
        d = min(maxval, death(i))
        b = birth(i)
        push!(B, (min(m*(d-b), b), d-b))
    end
    # pad zeros
    while length(B) < l+p
        push!(B, (zero(T), zero(T)))
    end
    n = length(B)

    if l == 0 && p == 0
        0
    elseif l == 0
        sort!(map(sum, B), rev=true)[1:p] |> sum
    elseif p == 0
        sort!(map(last,B), rev=true)[1:l] |> sum
    else
        res = 0
        L = hcat(combinations(1:n,l)...)'
        for i in 1:length(L)
            ii = L[i,:]
            e0 = sum(last, B[ii])
            K = hcat(combinations(setdiff(1:n, ii),p)...)'
            for j = 1:length(K)
                jj = K[j,:]
                e1 = sum(map(first, B[jj]) .+ map(last, B[jj]))
                res = max(res, e0 + e1)
            end
        end
        res
    end
end

"""
    tropic(pd)

Calculate a tropical coordinates for a persistent diagram `pd`

Paper: "Tropical Sufficient Statistics for Persistent Homology" by A.Monod et al.
Ref: https://arxiv.org/abs/1709.02647
"""
function tropic(pds::AbstractVector{PD}) where {T<:Real, PD<:PersistenceDiagram{T}}
    d,n,m,mv = ComputationalHomology.tropicstats(pds)
    [tropicp(pd,m,i,j,maxval=mv) for (i,j) in ComputationalHomology.Orbits, pd in pds ]
end
tropic(pd::PersistenceDiagram{T}) where {T<:Real} = tropic([pd]) |> vec
tropic(pdd::Dict{Int, PD}) where {T<:Real, PD<:PersistenceDiagram{T}} = tropic(flatten(pdd))
