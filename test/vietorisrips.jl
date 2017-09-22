@testset "Vietoris–Rips Complex" begin
    # Initialize dataset
    srand(9236493643764)
    N = 10
    X = rand(3,N)

    cplx, w = vietorisrips(X, 0.1)
    @test size(cplx, 0) == N
    @test size(cplx, 1) == 0

    cplx, w = vietorisrips(X, 0.3)
    @test size(cplx, 0) == N
    @test size(cplx, 1) == 4

    cplx, w = vietorisrips(X, 0.4, true)
    @test size(cplx, 0) == N
    @test size(cplx, 1) == 8
    @test size(cplx, 2) == 1
    @test w[0][1] == 0.
    @test w[1][1] ≈ 0.24932439624731
    @test w[2][1] ≈ 0.36949278019681

    cplx2, w2 = vietorisrips(X, 0.4, true, method=:inductive)
    @test size(cplx2, 0) == N
    @test size(cplx2, 1) == 8
    @test size(cplx2, 2) == 1
    @test w2[0][1] == 0.
    @test w2[1][1] ≈ 0.24932439624731
    @test w2[2][1] ≈ 0.36949278019681

    @test_throws ArgumentError vietorisrips(X, 0.4, true, method=:a)
end
