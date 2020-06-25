@testset "Boundaries" begin

    # boundary of a point
    ch = boundary(Simplex(1))
    @test iszero(ch)
    @test dim(ch) == -1

    # different order of faces
    @test Simplex(4,1) == Simplex(1,4)
    ch1 = boundary(Simplex(4,1))
    ch2 = boundary(Simplex(1,4))
    @test ch1 != ch2

    cplx = SimplicialComplex(Simplex(1,2,3),
                            Simplex(2,4),
                            Simplex(3,4),
                            Simplex(5,4),
                            Simplex(1,2,3,5))

    # ch = boundary(Simplex(1,2,3))
    # showchain(cplx, ch)

    # simplex boundary
    ch = Chain(1, Int) + (Simplex(1,2),1) + (Simplex(2,3),1) + (Simplex(1,3),-1)
    @testset "Complex boundary chain" for (a,b) in zip(boundary(Simplex(1,2,3)), ch)
        @test a == b
    end

    # ∂∘∂(c) == 0
    @test iszero( boundary(cplx, boundary(Simplex(1,2,3))) )

    # boundary of a path
    cplx = SimplicialComplex(Simplex(1,2), Simplex(2,3), Simplex(3,4))
    path = Chain(1, Dict(hash(c) => 1 for c in cells(cplx,1)))
    ch = Chain(0, Int) + (Simplex(1), -1) + (Simplex(4), 1)
    @testset "Path boundary" for (a,b) in zip(boundary(cplx, path), ch)
        @test a == b
    end

    # boundary of a loop
    push!(cplx, Simplex(4,1))
    loop = Chain(1, Dict(hash(c) => 1 for c in cells(cplx,1)))
    @test iszero(simplify(boundary(cplx, loop)))

    # coboundary
    ch = Chain(2, Int) + (Simplex(1,3,5),1) + (Simplex(1,2,3),-1)
    @testset "Complex coboundary chain" for (a,b) in zip(coboundary(cplx, Simplex(1,3)), ch)
        @test a == b
    end

    # δ∘δ(c) == 0
    @test iszero( coboundary(cplx, coboundary(cplx, Simplex(1,3))) )

    # cw-cell boundary
    a = Cell()
    b = Cell()
    c = Cell(1, a, b)
    ch = Chain(1, Int) + (a,1) + (b,-1)
    @testset "Complex boundary chain" for (a,b) in zip(boundary(c), ch)
        @test a == b
    end

    # cube boundary
    c = Cube([1,1]) # 2-cube
    ch = Chain(1, Int) + (Cube([0,1]), -1) + (Cube([1,0]), 1)
    @testset "Complex boundary chain" for (a,b) in zip(boundary(c), ch)
        @test a == b
    end

end