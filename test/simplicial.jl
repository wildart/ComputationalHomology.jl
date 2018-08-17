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
    @test size(cplx,0) == 4
    @test size(cplx,1) == 4
    @test size(cplx,2) == 1
    @test size(cplx,3) == 0
    @test cplx[Simplex(4), 0] == 4
    @test cplx[Simplex(10), 0] > size(cplx,0)
    @test cplx[Simplex(10), 0] == 5

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

    splx = push!(cplx, Simplex(1))
    @test splx[1][:index] == 1
    @test dim(splx[1]) == 0
    @test sum(size(cplx)) == 1
    @test cplx.cells[0][1][:values] == [1]

    splx = push!(cplx, Simplex(1,2))
    @test splx[1][:index] == 1
    @test dim(splx[1]) == 1
    @test sum(size(cplx)) == 2
    @test cplx.cells[1][1][:values] == [1,2]

    cplx = SimplicialComplex(Int)
    splxs = push!(cplx, Simplex(1,2), recursive=true)
    @test splxs[1][:index] == 1
    @test splxs[1][:values] == [1,2]
    @test splxs[2][:index] == 1
    @test splxs[2][:values] == [1]
    @test splxs[3][:index] == 2
    @test splxs[3][:values] == [2]
    @test sum(size(cplx)) == 3
    @test size(cplx, 0) == 2
    @test size(cplx, 1) == 1

    splxs = push!(cplx, Simplex(1,3), recursive=true)
    @test splxs[1][:index] == 2
    @test splxs[1][:values] == [1,3]
    @test splxs[2][:index] == 3
    @test splxs[2][:values] == [3]
    @test sum(size(cplx)) == 5
    @test size(cplx, 0) == 3
    @test size(cplx, 1) == 2

    io = IOBuffer()
    write(io, cplx)
    seekstart(io)
    cplx2 = read(io, SimplicialComplex{Int})
    @test size(cplx2, 0) == 3
    @test size(cplx2, 1) == 2

    ch = boundary(cplx2, 1, 0, Int)
    @test dim(ch) == -1
end
