# Persistent landscape

birth(i::Interval) = i[1] - i[2]
death(i::Interval) = i[1] + i[2]
Base.diag(i::Interval) = let c = death(i)/2.0; c => c end

rescaledrank(i::Interval) = (i[1] + i[2])/2.0 => (i[2] - i[1])/2.0
rescaledrank(bar::Vector{Interval}) = map(rescaledrank, sort(bar))
rescaledrank(bars::Dict{Int, Vector{Interval}}) = Dict(d => rescaledrank(bar) for (d,bar) in bars)

"""Construction of a persistent landscape from a barcode

Algorithm taken form "A Persistence Landscapes Toolbox for Topological Statistics" by  Bubenik & Dlotko, p. 22
"""
function landscape(bar::Vector{Interval})
    A = sort(bar, lt=(a,b) -> (a.first < b.first) ? true : (a.first > b.first ? false : (a.second > b.second)))
    k = 1
    L = Dict{Int, Vector{Interval}}()
    while length(A) > 0
        L[k] = Interval[]
        b, d  = A[1]
        deleteat!(A, 1)
        p = 1
        if (b, d) == (-Inf, Inf)
            push!(L[k],  -Inf => Inf)
        else
            if d == Inf
                push!(L[k], -Inf => 0.0)
                push!(L[k], b => 0.0)
                push!(L[k], Inf=>Inf)
            else
                if b == -Inf
                    push!(L[k], -Inf => Inf)
                else
                    push!(L[k], -Inf => 0.0)
                    push!(L[k], b => 0.0)
                    push!(L[k], rescaledrank(b=>d))
                end
            end
        end
        while !(L[k][end] == (Inf=>0.0) || L[k][end] == (Inf=>Inf))
            i = findfirst(x->last(x) > d, A[p:end])
            if i === nothing
                push!(L[k], d => 0.0)
                push!(L[k], Inf => 0.0)
            else
                b′, d′  = A[p+i-1]
                deleteat!(A, p+i-1)
                b′ > d && push!(L[k],  d => 0.0)
                if b′ >= d
                    push!(L[k], b′ => 0.0)
                else
                    push!(L[k], rescaledrank(b′=>d))
                    insert!(A, p+i-1, b′ => d)
                    p += 1
                end
                if d′ == Inf
                    push!(L[k], Inf=>Inf)
                else
                    push!(L[k], rescaledrank(b′=>d′))
                end
                b, d  = b′, d′
            end
        end
        k += 1
    end
    return L
end
