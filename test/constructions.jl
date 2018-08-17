@testset "Constructions" begin

    # Initialize dataset
    Random.seed!(9236493643764)
    N = 10
    X = rand(3,N)

    @testset "Vietoris–Rips Complex" begin
        cplx, w = vietorisrips(X, 0.1, false)
        @test size(cplx, 0) == N
        @test size(cplx, 1) == 0
        @test w === nothing

        cplx, w = vietorisrips(X, 0.3, false)
        @test size(cplx, 0) == N
        @test size(cplx, 1) == 4
        @test cplx[1,1] == Simplex(1, 10)
        @test w === nothing

        cplx, w = vietorisrips(X, 0.4)
        @test size(cplx, 0) == N
        @test size(cplx, 1) == 8
        @test cplx[1,1] == Simplex(1, 10)
        @test size(cplx, 2) == 1
        @test cplx[1,2] == Simplex(2, 4, 6)
        @test w[0][1] == 0.
        @test w[1][1] ≈ 0.24932439624731
        @test w[2][1] ≈ 0.36949278019681

        cplx2, w2 = vietorisrips(X, 0.4, expansion=:inductive)
        @test size(cplx2, 0) == N
        @test size(cplx2, 1) == 8
        @test size(cplx2, 2) == 1
        @test w2[0][1] == 0.
        @test w2[1][1] ≈ 0.24932439624731
        @test w2[2][1] ≈ 0.36949278019681

        @test_throws ArgumentError vietorisrips(X, 0.4, expansion=:a)
    end

    @testset "Witness Complex" begin
        l = 4

        cplx, w, L = witness(X, l, 0.005, false, firstpoint = 1)
        @test size(cplx, 0) == l
        @test size(cplx, 1) == 4
        @test L[cplx[1,1][:values]] == [1, 4]
        @test w === nothing

        cplx, w = witness(X, l, 0.1, false, firstpoint = 1)
        @test size(cplx, 0) == l
        @test size(cplx, 1) == 5
        @test L[cplx[5,1][:values]] == [4, 8]
        @test w === nothing

        cplx, w = witness(X, l, 0.5, firstpoint = 1)
        @test size(cplx, 0) == l
        @test size(cplx, 1) == 6
        @test L[cplx[1,1][:values]] == [1, 9]
        @test size(cplx, 2) == 4
        @test L[cplx[1,2][:values]] == [1, 9, 4]
        @test w[0][1] == 0.
        @test w[1][1] ≈ 0.22805075942971
        @test w[2][end] ≈ 0.0782128253900049

        cplx2, w2 = witness(X, l, 0.5, expansion=:inductive, firstpoint = 1)
        @test size(cplx2, 0) == l
        @test size(cplx2, 1) == 6
        @test L[cplx2[1,1][:values]] == [1, 9]
        @test size(cplx2, 2) == 4
        @test L[cplx2[1,2][:values]] == [1, 9, 4]
        @test w2[0][1] == 0.
        @test w2[1][1] ≈ 0.22805075942971
        @test w2[2][end] ≈ 0.0782128253900049

        @test_throws ArgumentError witness(X, l, 0.4, expansion=:a)

        L = ComputationalHomology.landmarks(X, l, method=:random)
        @test length(L) == l
    end
end
