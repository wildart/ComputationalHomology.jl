#=== Cube ===#

mutable struct Cube{P<:Vector} <: AbstractCell
    origin::P
    extent::P
    function Cube{P}(o::P, x::P) where {P}
        @assert length(x) <= length(o) "Too many extent elements"
        new(o, x)
    end
end
show(io::IO, c::Cube) = show(io, "Cube[$(c.origin) + $(c.extent)]")

# Public methods

dim(c::Cube) = length(findall(!iszero, c.extent))

==(a::Cube, b::Cube) = hash(a) == hash(b)

hash(c::Cube) = hash(vcat(c.origin, c.extent))

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
