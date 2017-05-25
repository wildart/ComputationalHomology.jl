# Homology: type and methods

abstract AbstractHomology{G}
grouptype{G}(::Type{AbstractHomology{G}}) = G
grouptype{H <: AbstractHomology}(::Type{H}) = supertype(H) |> grouptype
group{G}(h::AbstractHomology{G}, dim::Int) = throw(MethodError(group,(typeof(h),Int)))

"""
Homology group iterator for an abstract complex
"""
immutable Homology{C<:AbstractComplex, G} <: AbstractHomology{G}
    complex::C
end
homology{C <: AbstractComplex, G}(c::C, ::Type{G}) = Homology{C,G}(c)
Base.show(io::IO, h::Homology) = print(io, "Homology[$(h.complex)]")

"""Return homology group type: dimension, betti & torsion numbers."""
Base.eltype{C,G}(::Type{Homology{C,G}}) = Tuple{Int, Int, Int}

function group{C <: AbstractComplex, G}(h::Homology{C, G}, d::Int, torsionCoeffPrev::Int=0)
    M = boundary_matrix(G, h.complex, d+1)
    U, Uinv, V, Vinv, D = SNF(M)

    s = 0
    t = 0
    c = minimum(size(D))

    # Calculate rank B_{p−1} = rank D_p
    for i in 1:c
        if D[i,i] == one(G) || D[i,i] == -one(G)
            s += 1
        end
    end
    # Calculate rank Z_p = n_p − rank D_p
    for i in 1:c
        D[i,i] == zero(G) && break
        t += 1
    end

    # β_p = rank H p = rank Z_p - rank B_p = n_p − rank D_p - rank D_{p+1}
    bettiNum = size(U,1) - torsionCoeffPrev - t
    torsionCoeff = t - s

    return (bettiNum, torsionCoeff, (U, Uinv, V, Vinv, D, t))
end

#
# Iterator methods
#
Base.length{C, G}(h::Homology{C, G}) = dim(h.complex)

function Base.start{C, G}(h::Homology{C, G})
    Z = zeros(G,0,0)
    return (0, 0, (Z,Z,Z,Z,Z,0))
end

function Base.next{C, G}(h::Homology{C, G}, state)
    d = state[1]
    torsionCoeffPrev = state[2][end]

    bettiNum, torsionCoeff, snfstate = group(h, d, torsionCoeffPrev)

    return (d, bettiNum, torsionCoeff), (d+1, snfstate)
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
"""
withgenerators{H <: AbstractHomology}(h::H) = WithGenerators{H}(h)
Base.eltype{H <: AbstractHomology}(::Type{WithGenerators{H}}) = Tuple{eltype(H), Dict{Chain, grouptype(H)}}

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
    d = itm[1]
    bettiNum = itm[2]
    torsionCoeff  = itm[3]
    chains = Dict{Chain, GT}()
    U, Uinv, V, Vinv, D, t = nextState[2]
    cplx = g.homology.complex

    A = view(U, 1:size(U,1), t+1:size(U,2))
    B = view(Uinv, t+1:size(Uinv,1), 1:size(Uinv,2))
    C = if d == 0
        n = size(cplx, d)
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
        c = Chain(d, GT)
        for (i,v) in zip(I,V)
            push!(c, v=>i)
        end
        chains[c] = 0
    end

    # torsion generator
    for j in 0:torsionCoeff-1
        c = Chain(d, GT)
        I, V = findnz(U[:,s+j])
        for (i,v) in zip(I,V)
            push!(c, v=>i)
        end
        chains[c] = D[s+j, s+j]
    end

    return (itm, chains), nextState
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
    for ((d,b,t), gens) in g
        chains[d] = collect(keys(gens))
    end
    return chains
end
