@testset "Landscape" begin

    import ComputationalHomology: scale, value

    p = 1.0 => 2.0
    q = 0.0 => 3.0

    @test scale(q) == (1.5 => 1.5)
    @test value(p, q, 1) == 2.0

    bar = [p,q]
    @testset for (i,j) in zip(scale(bar), [1.5=>1.5, 1.5=>0.5])
        @test i == j
    end

    l = Landscape([p,q], [0.0=>1.0, 2.0=>Inf])
    @test length(l) == 2
    @test l[1][1] === p

    bar = [1=>5, 2=>8, 3=>4, 6=>7, 5=>9]
    L = Landscape(
        [ (-Inf=>0.0), (1=>0.0), (3.0=>2.0), (3.5=>1.5), (5.0=>3.0), (6.5=>1.5), (7.0=>2.0), (9=>0.0), (Inf=>0.0) ],
        [ (-Inf=>0.0), (2=>0.0), (3.5=>1.5), (5=>0.0), (6.5=>1.5), (8=>0.0), (Inf=>0.0) ],
        [ (-Inf=>0.0), (3=>0.0), (3.5=>0.5), (4=>0.0), (6=>0.0), (6.5=>0.5), (7=>0.0), (Inf=>0.0) ]
    )
    Lk = landscape(bar)
    @testset for (l1, l2) in zip(L, Lk)
        for (i1, i2) in zip(l1, l2)
            @test i1.first == i2.first
            @test i1.second == i2.second
        end
    end


    # reductions on landscapes
    l1 = Landscape([(-Inf=>0.0), (1.0=>0.0), (2.0=>1.0), (3.0=>0.0), (Inf=>0.0)])
    s = reduce(+, l1, 1.0)
    @testset for (p,q) in zip(l1[1], s[1])
        @test p.first == q.first
        @test p.second+1 == q.second
    end

    s = reduce(+, l1, l1)
    @test s[1][3] == (2 => 2)
    s = reduce(-, l1, l1)
    @test s[1][3] == (2 => 0)

    l2 = Landscape([(-Inf=>0.0), (2.0=>0.0), (3.0=>1.0), (4.0=>1.0), (5=>0.0), (Inf=>0.0)])
    s = reduce(+, l1, l2)
    @testset for (i, v) in zip(2:4, [(1 => 0), (2 => 1), (3 => 1)])
        @test s[1][i] == v
    end

    l3 = Landscape(l2[1], [(-Inf=>0.0), (1.0=>0.0), (2.0=>2.0), (3.0=>0.0), (Inf=>0.0)])
    s = reduce(+, l1, l3)
    @testset for (i, v) in zip(2:4, [(1 => 0), (2 => 1), (3 => 1)])
        @test s[1][i] == v
    end
    @test s[2][3] == (2 => 2)

    s = reduce(*, l3, l1)
    @testset for (i, v) in zip(2:4, [(1 => 0), (2 => 0), (3 => 0)])
        @test s[1][i] == v
    end
    @test s[2][3] == (2 => 0)

    m = mean([l1, l2, l3])
    @testset for (i, v) in zip(3:5, [(2 => 1/3), (3 => 2/3), (4 => 2/3)])
        @test m[1][i] == v
    end
    @test m[2][3] == (2 => 2/3)

end
