@testset "Boundaries" begin

    cplx = SimplicialComplex(Simplex(1,2,3),
                            Simplex(2,4),
                            Simplex(3,4),
                            Simplex(5,4),
                            Simplex(1,2,3,5))

    # ch = boundary(Simplex(1,2,3))
    # showchain(cplx, ch)

    # boundary of a point
    ch = boundary(Simplex(1))
    @test iszero(ch)
    @test dim(ch) == -1

    # simplex boundary
    ch = Chain(1, Int) + (Simplex(1,2),1) + (Simplex(2,3),1) + (Simplex(1,3),-1)
    @testset "Complex boundary chain" for (a,b) in zip(boundary(Simplex(1,2,3)), ch)
        @test a == b
    end

    # ∂∘∂(c) == 0
    @test iszero( boundary(cplx, boundary(Simplex(1,2,3))) )

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

end