@testset "Simplicial Complex" begin
    cplx = SimplicialComplex(Simplex(1,2,3), Simplex(2,4), Simplex(4))

    @test celltype(cplx) == Simplex{Int}
    @test length(cells(cplx)) == 3
    cs = cells(cplx,2)
    @test !isnull(cs)
    @test length(get(cs)) == 1
    cs = cells(cplx,3)
    @test isnull(cs)
    @test size(cplx) == (4,4,1)
    @test size(cplx,0) == 4
    @test size(cplx,1) == 4
    @test size(cplx,2) == 1
    @test size(cplx,3) == 0
    @test cplx[Simplex(4), 0] == 4
    @test cplx[Simplex(10), 0] > size(cplx,0)
    @test cplx[Simplex(10), 0] == 5

    cplx = SimplicialComplex(Simplex(1,2,3),
                       Simplex(2,4),
                       Simplex(3,4),
                       Simplex(5,4),
                       Simplex(6))

    qch = Chain(2, [1], [1])
    resch = Chain(1, [1, -1, 1], [3, 2, 1])
    @testset "Complex boundary chain" for (a,b) in zip(boundary(cplx, qch), resch)
        @test a == b
    end

    qch = Chain(1, [1], [1])
    resch = Chain(2, [-1, -1], [1, 2])
    @testset "Complex coboundary chain" for (a,b) in zip(coboundary(cplx, qch), resch)
        @test a == b
    end

    @test convert(Matrix, boundary_matrix(Int, cplx, 2)) == [1 -1 1 0 0 0]'

    cplx = SimplicialComplex(Int)
    s = addsimplex(cplx, Simplex(1))
    @test s == (1,0,1)
    @test cplx.cells[0][1][:values] == [1]
    s = addsimplex(cplx, Simplex(1,2))
    @test s == (2,1,1)
    @test cplx.cells[1][1][:values] == [1,2]

    splx = Simplex(1,2)
    cplx = SimplicialComplex(Int)
    s = addsimplex!(cplx, Simplex(1,2))
    @test s[1] == (1,0,1)
    @test s[2] == (2,0,2)
    @test s[3] == (3,1,1)
    s = addsimplex!(cplx, Simplex(1,3))
    @test s[1] == (4,0,3)
    @test s[2] == (5,1,2)
end
