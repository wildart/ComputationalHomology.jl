#=== Cube ===#

struct Cube{P} <: AbstractCell
    origin::AbstractVector{P}
    extent::AbstractVector{P}
    function Cube(o::AbstractVector{P}, x::AbstractVector{P}) where{P}
        @assert length(x) <= length(o) "Too many extent elements"
        new{P}(o, x)
    end
end
Cube(x::AbstractVector{P}) where {P} = Cube(zeros(P, length(x)), x)
Cube(::Type{P}) where {P} = Cube(P[], P[])
show(io::IO, c::Cube) = show(io, "Cube[$(c.origin) + $(c.extent)]")

# Public methods

dim(c::Cube) = length(findall(!iszero, c.origin - c.extent))

==(a::Cube, b::Cube) = hash(a) == hash(b)

hash(c::Cube) = hash(vcat(c.origin, c.extent))

function faces(c::Cube{P}, ε=zero(P)) where {P}
    d = dim(c)
    d == 0 && return Cube[]
    kₙ = findall(!iszero, c.origin - c.extent) # on-degenerate intervals
    (Cube(copy(c.origin), [i == j ? ε : c.extent[j] for j in 1:length(c.origin)]) for i in kₙ)
end

function boundary(::Type{R}, c::Cube) where {R}
    d = dim(c)
    ch = Chain(d-1, R)
    d == 0 && return ch
    o = one(R)
    for (i,(λ⁰, λ¹)) in enumerate(zip(faces(c,zero(R)), faces(c, one(R))))
        # println("λ⁰:", λ⁰," λ¹:", λ¹)
        push!(ch, hash(λ⁰)=>( isodd(i) ? -o : o))
        push!(ch, hash(λ¹)=>( isodd(i) ? o : -o))
    end
    return simplify(ch)
end

# Misc. methods

function volume(c::Cube)
    sizes = similar(c.origin,0)
    for (i,ex) in enumerate(c.extent)
        if ex != 0.
            expoint = copy(c.origin)
            expoint[i] += ex
            l = norm(expoint - c.origin)
            push!(sizes, l)
        end
    end
    return prod(sizes)
end
