@testset "Landscape" begin

    import ComputationalHomology: birth, death, rescaledrank

    p = 1.0 => 2.0
    q = 0.0 => 3.0

    @test birth(p) == -1.
    @test death(p) == 3.
    @test diag(p)  == (1.5 => 1.5)
    @test rescaledrank(q) == (1.5 => 1.5)

    bar = Interval[p,q]
    @testset for (i,j) in zip(rescaledrank(bar), Interval[1.5=>1.5, 1.5=>0.5])
        @test i == j
    end

    bars = Dict( 0 => Interval[p,q], 1 => Interval[0.0=>1.0, 2.0=>Inf] )
    @testset for (i,j) in zip(rescaledrank(bars)[1], Interval[0.5=>0.5, Inf=>Inf])
        @test i == j
    end

    bar = Interval[b=>d for (b,d) in [(1,5),(2,8),(3,4),(5,9),(6,7)]]
    L = [
        [ (-Inf,0.0), (1,0.0), (3.0,2.0), (3.5,1.5), (5.0,3.0), (6.5,1.5), (7.0,2.0), (9,0.0), (Inf,0.0)],
        [ (-Inf,0.0), (2,0.0), (3.5,1.5), (5,0.0), (6.5,1.5), (8,0.0), (Inf,0.0) ],
        [ (-Inf,0.0), (3,0.0), (3.5,0.5), (4,0.0), (6,0.0), (6.5,0.5), (7,0.0), (Inf,0.0)]
    ]
    Lk = landscape(bar)
    @testset for (k,l) in enumerate(L)
        for ((n1,r1),(n2,r2)) in zip(Lk[k], l)
            @test n1 == n2
            @test r1 == r2
        end
    end

end
