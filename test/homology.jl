@testset "Homology" begin
    cplx = SimplicialComplex(Simplex(1,2,3),
                        Simplex(2,4),
                        Simplex(3,4),
                        Simplex(5,4),
                        Simplex(6))

    h = homology(cplx)
    @test eltype(h) == Tuple{Int64,Int64,Int64}
    @test grouptype(supertype(typeof(h))) == Int
    @test grouptype(typeof(h)) == Int
    @test length(h) == 3

    b, t, _ = group(h,0);
    @test b == 2
    @test t == 0

    st = start(h)
    @test st[1] == 0

    itm, st = next(h, st)
    @test itm[1] == 0
    @test itm[2] == 2
    @test st[1] == 1

    itm, st = next(h, st)
    @test !done(h,st)
    @test itm[1] == 1
    @test itm[2] == 1
    @test st[1] == 2

    itm, st = next(h, st)
    @test done(h,st)
    @test itm[1] == 2
    @test itm[2] == 0
    @test st[1] == 3

    @test ComputationalHomology.betti(h) == [2,1,0]
    @test ComputationalHomology.euler(h) == 1

    g = withgenerators(h)
    @test eltype(g) == Tuple{Int64,Int64,Int64,Dict{Chain,Int64}}
    @test length(g) == 3

    st = start(g)
    @test st[1] == 0

    itm, st = next(g, st)
    @test itm[1] == 0
    @test itm[2] == 2
    @testset for (ch, d) in itm[4]
        @test d == 0
        @test ch[1][1] == 1
        @test ch[1][2] in [6,5]
    end

    itm, st = next(g, st)
    @test itm[1] == 1
    @test itm[2] == 1
    @test first(itm[4])[2] == 0
    @testset for (a,b) in zip( simplify(first(itm[4])[1]), Chain(1, Int) + (-1, 4) + (1,3) + (1,5))
        @test a == b
    end

    itm, st = next(g, st)
    @test itm[1] == 2
    @test itm[2] == 0
    @test length(itm[4]) == 0

    @test length(generators(g)) == 3
end
