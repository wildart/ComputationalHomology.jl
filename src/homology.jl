# Homology: type and methods

abstract AbstractHomology{G}
grouptype{G}(::Type{AbstractHomology{G}}) = G
grouptype{H <: AbstractHomology}(::Type{H}) = supertype(H) |> grouptype
group{G}(h::AbstractHomology{G}, dim::Int; kw...) = throw(MethodError(group,(typeof(h),Int)))

"""
Homology group iterator for an abstract complex
"""
immutable Homology{C<:AbstractComplex, G} <: AbstractHomology{G}
    complex::C
end
homology{C <: AbstractComplex, G}(c::C, ::Type{G}) = Homology{C,G}(c)
homology{C <: AbstractComplex}(c::C) = homology(c, Int)

Base.show(io::IO, h::Homology) = print(io, "Homology[$(h.complex)]")

"""Return homology group type: dimension, Betti & torsion numbers."""
Base.eltype{C,G}(::Type{Homology{C,G}}) = Tuple{Int, Int, Int}

#
# Interface methods
#

function group{C <: AbstractComplex, G}(h::Homology{C, G}, p::Int; Dₚ::Int=0)
    M = boundary_matrix(G, h.complex, p+1)
    U, Uinv, V, Vinv, D = SNF(M)

    nₚ, nₚ₊₁ = size(D)
    Dₚ₊₁ = trivial = 0

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

    return (βₚ, τₚ, (U, Uinv, V, Vinv, D, Dₚ₊₁))
end

#
# Iterator methods
#

Base.length{C, G}(h::Homology{C, G}) = dim(h.complex)+1

function Base.start{C, G}(h::Homology{C, G})
    Z = zeros(G,0,0)
    return (0, (Z,Z,Z,Z,Z,0))
end

function Base.next{C, G}(h::Homology{C, G}, state)
    p = state[1]

    βₚ, τₚ, snfstate = group(h, p, Dₚ = state[2][end])

    return (p, βₚ, τₚ), (p+1, snfstate)
end

function Base.done{C, G}(h::Homology{C, G}, state)
     return dim(h.complex) < state[1]
end


"""
Homology generator iterator
"""
immutable WithGenerators{H <: AbstractHomology}
    homology::H
end

"""
    withgenerators(hom)

An iterator that yields homology generators from a homology iterator `hom`.

Returns homology group parameters from iterator `hom` and generators as pair of `Chain` **x** and coefficient **k** that compose boundary **kx**.
When **k** is zero, then **x** is a boundary without any coefficient.
"""
withgenerators{H <: AbstractHomology}(h::H) = WithGenerators{H}(h)
Base.eltype{H <: AbstractHomology}(::Type{WithGenerators{H}}) = Tuple{Int, Int, Int, Dict{Chain, grouptype(H)}}

#
# Iterator methods
#
Base.length(g::WithGenerators) = length(g.homology)
Base.start(g::WithGenerators) = start(g.homology)
Base.done(g::WithGenerators, state) = done(g.homology, state)

function Base.next(g::WithGenerators, state)

    # calculate betti, torsion & SNF
    itm, nextState = next(g.homology, state)

    GT = grouptype(typeof(g.homology))
    p, βₚ, τₚ = itm
    chains = Dict{Chain, GT}()
    U, Uinv, V, Vinv, D, t = nextState[2]
    cplx = g.homology.complex
    trivial = t - τₚ

    A = view(U, 1:size(U,1), t+1:size(U,2))
    B = view(Uinv, t+1:size(Uinv,1), 1:size(Uinv,2))
    C = if p == 0
        n = size(cplx, p)
        speye(GT, n, n)
    else
        VinvPrev = state[2][4]
        tPrev = state[2][6]
        view(VinvPrev, 1:size(VinvPrev,1), tPrev+1:size(VinvPrev,2))
    end

    P = A * (B * C)

    gU, gUinv, gV, gVinv, gD = SNF(P)
    G = gU * gD

    # betti generator
    for j in 1:size(G,2)
        I, V = findnz(G[:,j])
        length(I) == 0 && continue
        ch = Chain(p, GT)
        for (i,v) in zip(I,V)
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

    return (p, βₚ, τₚ, chains), nextState
end

#
# Auxiliary  methods
#

"""Calculate Betti number"""
betti(g::AbstractHomology) = [gst[2] for gst in g]

"""Calculate Euler characteristic"""
euler(g::AbstractHomology) = [isodd(i) ? gst[2] : -gst[2] for (i,gst) in enumerate(g)] |> sum

"""Return linearly independent generators"""
function generators{H <: AbstractHomology}(g::WithGenerators{H})
    chains = Dict{Int,Vector{Chain}}()
    for (d,b,t,gens) in g
        chains[d] = collect(keys(gens))
    end
    return chains
end
