@testset "Persistent Image" begin

    intrs = intervals(0, 0=>Inf, 0=>1, 1=>1, 1=>2)

    pimg = PersistentImage(intrs...)
    @test length(pimg.intervals) == length(intrs)
    @test pimg(0,0) ≈ 0.62 atol=1-e2
    @test length(vec(pimg, 5, 6)) == 30

end

