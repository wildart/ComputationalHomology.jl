# Persistent image

struct PersistentImage
    intervals::Vector{Interval}
    weight::Function
    σ::Real
end
PersistentImage(ints::Interval...; weight=(b,d)->atan((d-b)/2.0), σ::Real = 2.0) =
    PersistentImage([ints...], weight, σ)

Base.show(io::IO, ρ::PersistentImage) = print(io, "ρ[$(length(ρ.intervals))]")

function (pimg::PersistentImage)(x, y)
    return sum(pimg.weight(i.b, i.d)*exp(-((i.b-x)^2 + (i.d-y)^2)/2*pimg.σ^2) for i in pimg.intervals)
end

function Base.vec(pimg::PersistentImage, B::Int, D::Int)
    intr = filter(!isinf∘birth, pimg.intervals)
    bmin, bmax = extrema(map(i->i.b, intr))
    dmin, dmax = extrema(map(i->i.d, intr))
    return [pimg(x,y) for x in range(bmin, bmax, length=B), y in range(dmin, dmax, length=D)] |> vec
end

