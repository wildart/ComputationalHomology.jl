# Persistent landscape

scale(p::Pair) = (p.first + p.second)/2 => (p.second - p.first)/2
scale(ps::Vector{<:Pair}) = map(scale, sort(ps))

struct Landscape
    levels::Vector{Vector{<:Pair}}
end
Landscape(ints::Vector{<:Pair}...) = Landscape([ints...])
getindex(ls::Landscape, i::Integer) = ls.levels[i]
length(ls::Landscape) = length(ls.levels)
show(io::IO, ls::Landscape) = print(io, "λ[k = $(length(ls.levels))]")
function Base.dump(io::IOContext, ls::Landscape; kwargs...)
    print(io, "Landscape:")
    for (i, l) in enumerate(ls.levels)
        print(io, "\n\tk = $i: ")
        for itr in l
            print(io, "$itr ")
        end
    end
end
function iterate(ls::Landscape, state=nothing)
    state === nothing && return iterate(ls, 1)
    state > length(ls) && return nothing
    return ls[state], state+1
end

function value(i1::Pair, i2::Pair, x)
    a = (i2.second - i1.second)/(i2.first - i1.first)
    b = i1.second - (a == 0 ? 0 : a*i1.first)
    return a*x+b
end

"""Construction of a persistent landscape from a barcode

Algorithm is adopted form "A Persistence Landscapes Toolbox for Topological Statistics" by  Bubenik & Dlotko, p. 22
"""
function landscape(bar::Vector{<:Pair})
    A = sort(bar)
    k = 0
    L = Vector{Vector{Pair}}()
    while length(A) > 0
        k += 1
        push!(L, Vector{Pair}())
        b, d  = A[1]
        deleteat!(A, 1)
        p = 1
        if (b, d) == (-Inf, Inf)
            push!(L[k],  -Inf => Inf)
        else
            if d == Inf
                push!(L[k], -Inf  => 0.0)
                push!(L[k], b => 0.0)
                push!(L[k], Inf  => Inf)
            else
                if b == -Inf
                    push!(L[k], -Inf => Inf)
                else
                    push!(L[k], -Inf => 0.0)
                    push!(L[k], b => 0.0)
                    push!(L[k], scale(b => d))
                end
            end
        end
        while !(L[k][end] == (Inf => 0.0) || L[k][end] == (Inf => Inf))
            i = findfirst(x->x.second > d, A[p:end])
            if i === nothing
                push!(L[k], d => 0.0)
                push!(L[k], Inf => 0.0)
            else
                b′, d′  = A[p+i-1]
                deleteat!(A, p+i-1)
                b′ > d && push!(L[k], d => 0.0)
                if b′ >= d
                    push!(L[k], b′ => 0.0)
                else
                    push!(L[k], scale(b′ => d))
                    insert!(A, p+i-1, b′ => d)
                    p += 1
                end
                if d′ == Inf
                    push!(L[k], Inf => Inf)
                else
                    push!(L[k], scale(b′ => d′))
                end
                b, d  = b′, d′
            end
        end
    end
    return Landscape(L)
end

"""
    reduce(op, l1::Landscape, x::Number)

Perform a reduction `op` between landscape and number.
"""
function Base.reduce(op, l::Landscape, x::Number)
    n = length(l)
    L = Vector{Vector{Pair}}(undef, n)
    for i in 1:n
        λᵢ = Pair[]
        for j in 1:length(l[i])
            push!(λᵢ, l[i][j].first => op(l[i][j].second, x))
        end
        L[i] = λᵢ
    end
    return Landscape(L);
end

"""
    reduce(op, l1::Landscape, l2::Landscape)

Perform a reduction `op` of two persistence landscapes.
"""
function Base.reduce(op, l1::Landscape, l2::Landscape)
    n = min(length(l1), length(l2))
    m = max(length(l1), length(l2))
    L = Vector{Vector{Pair}}(undef, m)
    for i in 1:n
        λᵢ = Pair[]
        p = q = 1
        while p+1 <= length(l1[i]) && q+1 <= length(l2[i])
            if l1[i][p].first < l2[i][q].first
                v = value(l2[i][q-1], l2[i][q], l1[i][p].first)
                push!(λᵢ, l1[i][p].first => op(l1[i][p].second, v))
                @debug "L" p=p q=q λᵢ=λᵢ[end] val=v
                p+=1
            elseif l1[i][p].first > l2[i][q].first
                v = value(l1[i][p-1], l1[i][p], l2[i][q].first)
                push!(λᵢ, l2[i][q].first => op(v, l2[i][q].second))
                @debug "G" p=p q=q λᵢ=λᵢ[end] val=v
                q+=1
            else # l1[i][p].first == l2[i][q].first
                push!(λᵢ, l2[i][q].first => op(l1[i][p].second, l2[i][q].second))
                @debug "E" p=p q=q λᵢ=λᵢ[end]
                p+=1
                q+=1
            end
        end
        while p+1 < length(l1[i]) && q+1 >= length(l2[i])
            push!(λᵢ, l1[i][p].first => op(l1[i][p].second, 0))
            @debug "L1" p=p q=q λᵢ=λᵢ[end]
            p+=1
        end
        while p+1 >= length(l1[i]) && q+1 < length(l2[i])
            push!(λᵢ, l2[i][q].first => op(0, l2[i][q].second))
            @debug "L2" p=p q=q λᵢ=λᵢ[end]
            q+=1
        end
        push!(λᵢ, Inf => 0.0)
        L[i] = λᵢ
    end
    if length(l1) > n
        for i in n+1:m
            λᵢ = copy(l1[i])
            for j in 1:length(λᵢ)
                λᵢ[j] = l1[i][j].first => op(l1[i][j].second, 0)
            end
            L[i] = λᵢ
            @debug "L1N:" i=i n=n m=m λᵢ=λᵢ
        end
    end
    if length(l2) > n
        for i in n+1:m
            λᵢ = copy(l2[i])
            for j in 1:length(λᵢ)
                λᵢ[j] = l2[i][j].first => op(0, l2[i][j].second)
            end
            L[i] = λᵢ
            @debug "L2N:" i=i n=n m=m λᵢ=λᵢ
        end
    end
    return Landscape(L)
end

"""
    mean(ls::Vector{Landscape})

Calculate avarage landscape.
"""
function mean(ls::Vector{Landscape})
    n = length(ls)
    m = foldl((l,v) -> reduce(+, l, v), ls)
    return reduce(/, m, n)
end
