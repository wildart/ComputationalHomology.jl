@testset "Simplicial Complex" begin
    sc = SimplicialComplex(Simplex(1,2,3), Simplex(2,4), Simplex(4))

    @test celltype(sc) == Simplex{Int}
    @test length(cells(sc)) == 3
    cs = cells(sc,2)
    @test !isnull(cs)
    @test length(get(cs)) == 1
    cs = cells(sc,3)
    @test isnull(cs)
    @test size(sc) == (4,4,1)
    @test size(sc,0) == 4
    @test size(sc,1) == 4
    @test size(sc,2) == 1
    @test size(sc,3) == 0
    @test sc[Simplex(4), 0] == 4
    @test sc[Simplex(10), 0] > size(sc,0)
    @test sc[Simplex(10), 0] == 5

    sc = SimplicialComplex(Simplex(1,2,3),
                       Simplex(2,4),
                       Simplex(3,4),
                       Simplex(5,4),
                       Simplex(6))

    qch = Chain(2, [1], [1])
    resch = Chain(1, [1, -1, 1], [3, 2, 1])
    @testset "Complex boundary chain" for (a,b) in zip(boundary(sc, qch), resch)
        @test a == b
    end

    qch = Chain(1, [1], [1])
    resch = Chain(2, [-1, -1], [1, 2])
    @testset "Complex coboundary chain" for (a,b) in zip(coboundary(sc, qch), resch)
        @test a == b
    end

    @test convert(Matrix, boundary_matrix(Int, sc, 2)) == [1 -1 1 0 0 0]'

    sc = SimplicialComplex(Int)
    s = addsimplex(sc, Simplex(1))
    @test s == (1,0,1)
    @test sc.cells[0][1][:values] == [1]
    s = addsimplex(sc, Simplex(1,2))
    @test s == (2,1,1)
    @test sc.cells[1][1][:values] == [1,2]

    sc = SimplicialComplex(Int)
    s = addsimplex!(sc, Simplex(1,2))
    @test s[1] == (1,0,1)
    @test s[2] == (2,0,2)
    @test s[3] == (3,1,1)
    s = addsimplex!(sc, Simplex(1,3))
    @test s[1] == (4,0,3)
    @test s[2] == (5,1,2)
end
