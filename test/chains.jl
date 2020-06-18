@testset "Chains" begin
    a = Chain(Int8)
    @test dim(a) == 0
    @test length(a) == 0
    @test keytype(a) == UInt
    @test valtype(a) == Int8

    a = Chain(Int, Int)
    @test keytype(a) == Int
    @test valtype(a) == Int

    b = Chain(Int[],Int[])
    @test keytype(b) == Int
    @test valtype(b) == Int
    push!(b, (1, 2))
    @test first(keys(b)) == 1
    @test b[1] == 2
    push!(b, 1, -1)
    @test b[1] == 1
    @test (b + (1, -1))[1] == 0
    @test b[1] == 1

    c = Chain(Int,Int)
    push!(c, 0=>1)
    @test length(c) == 1
    @test first(keys(c)) == 0
    @test first(values(c)) == 1

    d = Chain(Int,UInt8)
    d+=(1,0x01)
    @test d[1] == 1
    @test (d+(1,0x01))[1] == 2
    @test d[1] == 1

    a+=c
    a+=b
    a+=(2,1)
    @testset "loop thrugh chain" for i in keys(a)
        @test a[i] == 1
    end

    b+=(2,1)
    @test b[2] == 1

    append!(c,b)
    @testset "two chain compare" for (aa,cc) in zip(a,c)
        @test cc == aa
    end

    @testset "sum" for (i,v) in a+c
        @test v == 2
    end

    @test (2+b+2)[1] == 5
    @test (2*b*2)[1] == 4
    @test (1+2*b*3+4)[1] == 11

    @test iszero(simplify(a-c))
end
