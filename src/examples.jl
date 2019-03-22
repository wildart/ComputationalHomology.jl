
"""
    sphere(n)

A minimal triangulation of the `n`-dimensional sphere.
"""
function sphere(n::Integer)
    s = Simplex(collect(1:n+1))
    return SimplicialComplex(s)
end

"""
    torus()

A minimal triangulation of the torus.

This is a simplicial complex with 7 vertices, 21 edges and
14 faces. It is the unique triangulation of the torus with
7 vertices, and has been found by Möbius in 1861.
"""
function torus()
    cplx = SimplicialComplex{Simplex{Int}}()
    for s in [[1, 2, 3], [2, 3, 5], [2, 4, 5], [2, 4, 7], [1, 2, 6],
              [2, 6, 7], [3, 4, 6], [3, 5, 6], [3, 4, 7], [1, 3, 7],
              [1, 4, 5], [1, 4, 6], [5, 6, 7], [1, 5, 7]]
        addsimplices!(cplx, Simplex(s))
    end
    return cplx
end


"""
    kleinbottle()

A minimal triangulation of the Klein bottle, as presented for example
in Davide Cervone's thesis.
"""
function kleinbottle()
    cplx = SimplicialComplex{Simplex{Int}}()
    for s in [[3, 4, 8], [2, 3, 4], [2, 4, 6], [2, 6, 8], [2, 5, 8],
              [3, 5, 7], [2, 3, 7], [2, 7, 1], [2, 5, 1], [3, 5, 1],
              [4, 5, 8], [4, 5, 7], [4, 6, 7], [6, 7, 1], [3, 6, 1], [3, 6, 8]]
        addsimplices!(cplx, Simplex(s))
    end
    return cplx
end

"""
    projectiveplane()

A minimal triangulation of the real projective plane.
"""
function projectiveplane()
    cplx = SimplicialComplex{Simplex{Int}}()
    for s in [[1, 2, 3], [1, 3, 4], [1, 2, 6], [1, 5, 6], [1, 4, 5],
              [2, 3, 5], [2, 4, 5], [2, 4, 6], [3, 4, 6], [3, 5, 6]]
        addsimplices!(cplx, Simplex(s))
    end
    return cplx
end

"""
    randomcomplex(n, d)

A random `d`-dimensional simplicial complex on `n` vertices.
"""
function randomcomplex(n, d)
    cplx = SimplicialComplex(map(Simplex, 1:n))
    if d+1 >= n
        addsimplices!(cplx, Simplex(collect(1:n)))
    else
        pts = collect(1:n)
        simsset = Set{Set{Int}}()
        while length(simsset) < n
            Random.shuffle!(pts)
            push!(simsset, Set(pts[1:d+1]))
        end
        for s in simsset
            addsimplices!(cplx, Simplex(s...))
        end
    end
    return cplx
end

