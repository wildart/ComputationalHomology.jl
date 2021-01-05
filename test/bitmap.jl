@testset "Bitmap Complex" begin
    # empty complex
    cplx = BitmapComplex()
    @test size(cplx) == (0,)
    @test size(cplx,0) == 0
    @test length(cplx) == 0
    @test dim(cplx) == 0
    @test eltype(cplx) == Float64
    cplx = BitmapComplex(Int)
    @test eltype(cplx) == Int

    # 2D complex
    img = [
        7 6  5 ;
        8 20 4;
        1 4  6 ;
    ]
    cplx = BitmapComplex(img)
    @test size(cplx) == (16, 24, 9)
    @test size(cplx,0) == 16
    @test length(cplx) == 49
    @test dim(cplx) == 2
    @test eltype(cplx) == Int

    @test ComputationalHomology.topcubes(cplx) == [9, 11, 13, 23, 25, 27, 37, 39, 41]

    @test ComputationalHomology.cellindex(cplx, 1, 1) == 9

    cdc  = [ComputationalHomology.celldim(cplx, i) for i in 1:length(cplx)]
    cdt = vec([0 1 0 1 0 1 0 1 2 1 2 1 2 1 0 1 0 1 0 1 0 1 2 1 2 1 2 1 0 1 0 1 0 1 0 1 2 1 2 1 2 1 0 1 0 1 0 1 0])
    @test cdc == cdt

    @test faces(cplx, 9) == [2, 16, 10, 8]
    @test cofaces(cplx, 3) == [10, 2, 4]

    flt = filtration(cplx)
    @test eltype(flt) == BitmapComplex{Int}
    @test valtype(flt) == Float32
    @test minimum(flt) == 1
    @test maximum(flt) == 20

    ∂ = boundary(flt)
    @test length(∂) == length(complex(flt))
    @test count(!iszero, sparse(∂)) == 84

    @testset "PH Intervals " for (d, itrs) in diagram(flt)
        titrs = d == 0 ? diagram(1.0=>Inf) : diagram(8.0=>20.0)
        for itr in itrs
            @test itr ∈ titrs
        end
    end

    # Perseus format
    str = "2\n3\n3\n1\n4\n6\n8\n20\n4\n7\n6\n5"
    io = IOBuffer(str)
    cplx = read(io, BitmapComplex, Int)
    @test size(cplx) == (16, 24, 9)
    @test dim(cplx) == 2
    @test eltype(cplx) == Int

    # img = reshape([0,0,2,2,2,3,0,0,0,2,0,0,1,3,0,0,0,0,2,1,2,0,0,0,0,0,0,4,1,0,0,0,0,0,1,1,1,1,0,0,0,0,2,0,0,2,1,0,0,0,2,0,0,0,2,0,0,0,1,2,2,3,0,0],8,8)
    # flt2 = filtration(BitmapComplex(maximum(img) .- img))
    # diagram(flt2)
    # z = maximum(img) .- img
    # join(mapslices(v->join(v,"\t"), map(v->v == 4 ? "" : string(v), z), dims=2), "\n") |> println

end
