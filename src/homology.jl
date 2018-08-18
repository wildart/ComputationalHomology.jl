# Homology: type and methods

abstract type AbstractHomology end
grouptype(::Type{AbstractHomology}) = Nothing
grouptype(::Type{H}) where {H <: AbstractHomology} = supertype(H) |> grouptype
group(h::AbstractHomology, dim::Int; kw...) = throw(MethodError(group,(typeof(h),Int)))

"""
Homology group iterator for an abstract complex
"""
struct Homology{C<:AbstractComplex, G} <: AbstractHomology
    complex::C
end
homology(c::C, ::Type{G}) where {C <: AbstractComplex,G} = Homology{C,G}(c)
homology(c::C) where {C <: AbstractComplex} = homology(c, Int)

grouptype(::Type{Homology{C,G}}) where {C,G} = G

Base.show(io::IO, h::Homology) = print(io, "Homology[$(h.complex)]")

"""Return homology group type: dimension, Betti & torsion numbers."""
Base.eltype(::Type{Homology{C,G}}) where {C,G} = Tuple{Int, Int, Int}

#
# Interface methods
#

function group(h::Homology{C, G}, p::Int; Dₚ::Int=0) where {C <: AbstractComplex, G}
    cdim = dim(h.complex)
    @assert cdim >= p "Cannot define $p-th homology group for $cdim-dimensional complex"

    M = boundary_matrix(G, h.complex, p+1)
    F = smith(M)
    D = diagm(F)

    nₚ, nₚ₊₁ = size(D)
    Dₚ₊₁ = trivial = 0

    # empty boundary
    nₚ₊₁ == 0 && return (0, 0, (F.S, F.Sinv, F.T, F.Tinv, D, Dₚ₊₁))

    # Calculate rank B_p = rank D_{p+1}
    for i in 1:nₚ₊₁
        D[i,i] == zero(G) && break
        Dₚ₊₁ += 1
    end

    # rank of torsion-free subgroup
    for i in 1:min(nₚ, nₚ₊₁)
        if D[i,i] == one(G) || D[i,i] == -one(G)
           trivial += 1
        end
    end

    # β_p = rank H p = rank Z_p - rank B_p = n_p − rank D_p - rank D_{p+1}
    βₚ = nₚ - Dₚ - Dₚ₊₁
    τₚ = Dₚ₊₁ - trivial

    return (βₚ, τₚ, (F.S, F.Sinv, F.T, F.Tinv, D, Dₚ₊₁))
end

#
# Iterator methods
#

Base.length(h::Homology{C, G}) where {C,G} = dim(h.complex)+1
Base.eltype(h::Homology{C, G}) where {C,G} = Tuple{Int,Int,Int}

function Base.iterate(h::Homology{C, G}, state=nothing) where {C,G}
    if state === nothing
        Z = zeros(G,0,0)
        snfstate = (Z,Z,Z,Z,Z,0)
        return iterate(h, (0, snfstate))
    end

    p = state[1]
    p > dim(h.complex) && return nothing

    βₚ, τₚ, snfstate = group(h, p, Dₚ = state[2][end])
    return ((p, βₚ, τₚ), (p+1, snfstate))
end

"""
Homology generator iterator
"""
struct WithGenerators{H <: AbstractHomology}
    homology::H
end

"""
    withgenerators(hom)

An iterator that yields homology generators from a homology iterator `hom`.

Returns homology group parameters from iterator `hom` and generators as pair of `Chain` **x** and coefficient **k** that compose boundary **kx**.
When **k** is zero, then **x** is a boundary without any coefficient.
"""
withgenerators(h::H) where {H <: AbstractHomology} = WithGenerators{H}(h)
Base.eltype(::Type{WithGenerators{H}}) where {H <: AbstractHomology} = Tuple{Int, Int, Int, Dict{Chain, grouptype(H)}}

#
# Iterator methods
#
Base.length(g::WithGenerators) = length(g.homology)
Base.eltype(g::WithGenerators) = Tuple{Int,Int,Int,Dict{Chain, grouptype(typeof(g.homology))}}

function Base.iterate(g::WithGenerators, state=nothing)
    h = g.homology

    # calculate betti, torsion & SNF
    result = iterate(h, state)
    result === nothing && return nothing

    GT = grouptype(typeof(h))
    chains = Dict{Chain, GT}()

    p, βₚ, τₚ = result[1]
    U, Uinv, V, Vinv, D, t = result[2][2]
    cplx = g.homology.complex
    trivial = t - τₚ

    A = view(U, 1:size(U,1), t+1:size(U,2))
    B = view(Uinv, t+1:size(Uinv,1), 1:size(Uinv,2))
    C = if p == 0
        n = size(cplx, p)
        SparseMatrixCSC{GT}(I, n, n)
    else
        VinvPrev = state[2][4]
        tPrev = state[2][6]
        view(VinvPrev, 1:size(VinvPrev,1), tPrev+1:size(VinvPrev,2))
    end

    P = A * (B * C)
    F = smith(P)
    G = F.S * diagm(F)

    # betti generator
    for j in 1:size(G,2)
        ii, V = findnz(G[:,j])
        length(ii) == 0 && continue
        ch = Chain(p, GT)
        for (i,v) in zip(ii,V)
            push!(ch, v=>i)
        end
        chains[ch] = 0
    end

    # torsion generator
    for j in 1:τₚ
        ch = Chain(p, GT)
        I, V = findnz(D[:,trivial+j])
        for (i,v) in zip(I,V)
            push!(ch, v=>i)
        end
        chains[ch] = D[trivial+j, trivial+j]
    end

    return (p, βₚ, τₚ, chains), result[2]
end

#
# Auxiliary  methods
#

"""Calculate Betti number"""
betti(g::AbstractHomology) = [gst[2] for gst in g]

"""Calculate Euler characteristic"""
euler(g::AbstractHomology) = [isodd(i) ? gst[2] : -gst[2] for (i,gst) in enumerate(g)] |> sum

"""Return linearly independent generators"""
function generators(g::WithGenerators{H}) where {H <: AbstractHomology}
    chains = Dict{Int,Vector{Chain}}()
    for (d,b,t,gens) in g
        chains[d] = collect(keys(gens))
    end
    return chains
end
