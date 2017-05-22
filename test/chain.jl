@testset "Chain" begin
    a = Chain(0, Int)
    @test dim(a) == 0
    setdim!(a,1)
    @test dim(a) == 1

    b = Chain(Int[],Int[])
    c = Chain(Int)

    push!(c, (1,0))
    @test c[1] == (1,0)

    push!(b, 1=>1)
    @test b[1] == (1,1)

    a+=c
    a+=b
    @test a[1] == (1,0)
    @test a[2] == (1,1)

    @testset "loop thrugh chain" for (chc,chel) in a
        @test chc == 1
    end

    append!(c,b)
    @testset "two chain compare" for (aa,cc) in zip(a,c)
        @test cc == aa
    end

    @testset "simplify" for (chc,chel) in simplify(a+c)
        @test chc == 2
    end

    @test (2*b)[1] == (2,1)
    @test (b*2)[1] == (2,1)

    @test length(simplify(a-c)) == 0
end
