@testset "Distances" begin

    dgm1 = intervals(0, 0.5=>1.0)
    dgm2 = intervals(0, 0.5=>1.1)
    dgm3 = intervals(0, 0.5=>1.0, 0.6=>1.1, 0.3 => Inf)

    m = convert(Matrix, dgm1)
    @test size(m) == (2,1)
    m = convert(Matrix, dgm3)
    @test size(m) == (2,2)
    m = convert(Matrix, dgm3, skipinf=false)
    @test size(m) == (2,3)

    @test wasserstein(dgm1, dgm1) ≈ 0.0
    @test wasserstein(dgm1, dgm2) ≈ 0.1
    @test wasserstein(dgm2, dgm2) ≈ 0.0
    @test wasserstein(dgm3, dgm2) ≈ 0.2
    @test wasserstein(dgm2, dgm3) ≈ 0.2
    @test wasserstein(intervals(0, 1.0=>2.0), Interval{Float64}[]) ≈ √2/2
    @test wasserstein(intervals(0, 1.0=>2.0), intervals(0, 0.0=>Inf)) ≈ √2/2
    @test wasserstein(intervals(0, 0.0=>Inf), intervals(0, 1.0=>2.0)) ≈ √2/2

end
