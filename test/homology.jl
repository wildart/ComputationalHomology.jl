using PersistentHomology
using Base.Test

@testset "Homology" begin

    cplx = SimplicialComplex(Simplex(1,2,3),
                        Simplex(2,4),
                        Simplex(3,4),
                        Simplex(5,4),
                        Simplex(6))

    h = homology(cplx, Int)
    @test eltype(typeof(h)) == Tuple{Int64,Int64,Int64}
    @test grouptype(supertype(typeof(h))) == Int
    @test grouptype(typeof(h)) == Int
    @test length(h) == 2

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

    g = withgenerators(h)
    @test eltype(typeof(g)) == Tuple{Tuple{Int64,Int64,Int64},Dict{PersistentHomology.Chain,Int64}}
    @test length(g) == 2

    st = start(g)
    @test st[1] == 0

    itm, st = next(g, st)
    @test itm[1][1] == 0
    @test itm[1][2] == 2
    @testset for (ch, d) in itm[2]
        @test d == 0
        @test ch[1][1] == 1
        @test ch[1][2] in [6,5]
    end

    itm, st = next(g, st)
    @test itm[1][1] == 1
    @test itm[1][2] == 1
    @test first(itm[2])[2] == 0
    @testset for (a,b) in zip( simplify(first(itm[2])[1]), Chain(1, Int) + (-1, 4) + (1,3) + (1,5))
        @test a == b
    end

    itm, st = next(g, st)
    @test itm[1][1] == 2
    @test itm[1][2] == 0
    @test length(itm[2]) == 0

    @test length(generators(g)) == 3

end
