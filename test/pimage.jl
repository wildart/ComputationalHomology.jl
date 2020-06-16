@testset "Persistent Image" begin

    intrs = diagram(0, 0.0=>Inf, 0.0=>1.0, 1.0=>1.0, 1.0=>2.0)

    pimg = PersistentImage(intrs...)
    @test length(pimg.intervals) == length(intrs)
    @test pimg(0,0) ≈ 0.0627 atol=1e-3
    @test length(vec(pimg, 5, 6)) == 30

end

