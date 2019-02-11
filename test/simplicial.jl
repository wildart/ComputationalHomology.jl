@testset "Simplicial Complex" begin
    cplx = SimplicialComplex(Simplex(1,2,3), Simplex(2,4), Simplex(4))

    @test celltype(cplx) == Simplex{Int}
    @test length(cells(cplx)) == 3
    cs = cells(cplx,2)
    @test cs !== nothing
    @test length(cs) == 1
    cs = cells(cplx,3)
    @test cs === nothing
    @test size(cplx) == (4,4,1)
    @test size(cplx,-1) == 0
    @test size(cplx,0) == 4
    @test size(cplx,1) == 4
    @test size(cplx,2) == 1
    @test size(cplx,3) == 0
    @test cplx[Simplex(4)] == hash(Simplex(4))
    @test position(cplx, Simplex(4)) == 4
    @test cplx[Simplex(10)] === nothing
    @test position(cplx, Simplex(10)) === nothing

    @test length(simplices(cplx)) == 9
    scitr = simplices(cplx, 1)
    @test length(scitr) == 4
    for s in simplices(cplx, 2)
        @test s == Simplex(1,2,3)
    end

    cplx = SimplicialComplex(Simplex(1,2,3),
                       Simplex(2,4),
                       Simplex(3,4),
                       Simplex(5,4),
                       Simplex(6))

    # boundary
    qch = Chain(2, [1], [hash(Simplex(1,2,3))])
    resch = Chain(1, [-1, 1, 1], map(hash, [Simplex(1,2), Simplex(2,3), Simplex(1,3)]))
    @testset "Complex boundary chain" for (a,b) in zip(boundary(cplx, qch), resch)
        @test a == b
    end

    # ∂∘∂(c) == 0
    @test iszero(boundary(cplx, boundary(cplx, hash(Simplex(1,2,3)), 2)))

    @test convert(Matrix, boundary(cplx, 2, Int)) == [1 -1 1 0 0 0]'

    # coboundary
    push!(cplx, Simplex(1,3,5), recursive=true)
    push!(cplx, Simplex(1,2,3,5), recursive=true)

    qch = Chain(1, [1], [hash(Simplex(1,3))])
    resch = Chain(2, [-1, 1], map(hash,[Simplex(1,3,5), Simplex(1,2,3)]))
    @testset "Complex coboundary chain" for (a,b) in zip(coboundary(cplx, qch), resch)
        @test a == b
    end

    # δ∘δ(c) == 0
    @test iszero(coboundary(cplx, coboundary(cplx, hash(Simplex(1,3)), 1)))

    cplx = SimplicialComplex(Simplex{Int})

    splx = push!(cplx, Simplex(1))
    @test position(cplx, splx[1]) == 1
    @test dim(splx[1]) == 0
    @test sum(size(cplx)) == 1
    @test values(cells(cplx, 0)[1]) == [1]

    splx = push!(cplx, Simplex(1,2))
    @test splx[1] == Simplex(1,2)
    @test dim(splx[1]) == 1
    @test sum(size(cplx)) == 2
    @test sort!(values(cells(cplx, 1)[1])) == [1,2]

    cplx = SimplicialComplex(Simplex{Int})
    splxs = push!(cplx, Simplex(1, 2), recursive=true)
    @test splxs[1] == Simplex(1,2)
    @test sort!(values(splxs[1])) == [1, 2]
    @test splxs[2] == Simplex(2)
    @test values(splxs[2]) == [2]
    @test splxs[3] == Simplex(1)
    @test values(splxs[3]) == [1]
    @test sum(size(cplx)) == 3
    @test size(cplx, 0) == 2
    @test size(cplx, 1) == 1

    splxs = push!(cplx, Simplex(1,3), recursive=true)
    @test splxs[1] == Simplex(1,3)
    @test sort!(values(splxs[1])) == [1,3]
    @test splxs[2] == Simplex(3)
    @test values(splxs[2]) == [3]
    @test sum(size(cplx)) == 5
    @test size(cplx, 0) == 3
    @test size(cplx, 1) == 2

    @test Simplex(1,3) in cplx
    @test first(cells(cplx,1)) in cplx
    @test !(Simplex(2,3) in cplx)
    @test !(Simplex(2,3,3) in cplx)

    io = IOBuffer()
    write(io, cplx)
    seekstart(io)
    cplx2 = read(io, SimplicialComplex{Simplex{Int}})
    @test size(cplx2, 0) == 3
    @test size(cplx2, 1) == 2

    ch = boundary(cplx2, 1, 0, Int)
    @test dim(ch) == -1
end
