@testset "Constructions" begin

    # Initialize dataset
    Random.seed!(9236493643764)
    N = 10
    X = rand(3,N)

    # Gaudi example of Rips complex
    # http://gudhi.gforge.inria.fr/doc/latest/group__rips__complex.html
    Y = [ 1. 1.
          7. 0.
          4. 6.
          9. 6.
          0. 14.
          2. 19.
          9. 17.]

    # Gaudi example of Cech complex
    # http://gudhi.gforge.inria.fr/doc/latest/group__cech__complex.html
    Z = [ 1. 0.
          0. 1.
          2. 1.
          3. 2.
          0. 3.
          3+sqrt(3) 3.
          1. 4.
          3. 4.
          2. 4. + sqrt(3)
          0. 4.
         -0.5 2.]

    @testset "Čech Complex" begin
        r = 1.0

        cplx, w = čech(Z', r; maxoutdim = 2)
        @test sum(size(cplx)) == 30

        f = filtration(cplx, w)
        @testset "Test from Gudhi" for ((v,ls), d) in
            zip(simplices(f), [0. 0.5 0.559017 0.707107 0.99999 1.0])
            @test v ≈ d atol=1e-4
        end

        cplx, w = cech(Z', r, false)
        @test size(cplx, 0) == 11
        @test size(cplx, 1) == 16
        @test size(cplx, 2) == 3
        @test w === nothing
    end

    @testset "Vietoris-Rips Complex" begin

        @testset for expansion in [:inductive, :incremental]
            cplx, w = vietorisrips(Y', 11.; maxoutdim = 1, expansion=expansion)
            @test size(cplx, 0) == 7
            @test size(cplx, 1) == 11
            @test length(cplx) == 18

            f = filtration(cplx, w)
            @testset "Test from Gudhi" for ((v,ls), d) in
                zip(simplices(f),
                    [0. 5. 5.38516 5.83095 6.08276 6.32456 6.7082 7.28011 8.94427 9.43398 9.48683 11])
                @test v ≈ d atol=1e-4
            end

            cplx, w = vietorisrips(Y', 11.; expansion=expansion)
            @test length(cplx) == 24
            @test cplx.cells[3][1] == Simplex(1,2,3,4)
            @test w[3][1] ≈ 9.43398 atol=1e-4
            @test cplx.cells[1][8] == Simplex(4,7)
            @test w[1][8] == 11.0
        end

        cplx, w = vietorisrips(X, 0.1, false)
        @test size(cplx, 0) == N
        @test size(cplx, 1) == 0
        @test w === nothing

        cplx, w = vietorisrips(X, 0.3, false)
        @test size(cplx, 0) == N
        @test size(cplx, 1) == 4
        @test cells(cplx, 1)[1] == Simplex(1, 10)
        @test w === nothing

        cplx, w = vietorisrips(X, 0.4)
        @test size(cplx, 0) == N
        @test size(cplx, 1) == 8
        @test cells(cplx, 1)[1] == Simplex(1, 10)
        @test size(cplx, 2) == 1
        @test cells(cplx, 2)[1] == Simplex(2, 4, 6)
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

        L = ComputationalHomology.landmarks(X, l, method=:random)
        @test length(L) == l

        L = ComputationalHomology.landmarks(X, l, firstpoint = 1)
        @test L[1] == 1
        @test all(i->i ∈ L, [8, 1, 4, 9])

        L = ComputationalHomology.landmarks(X, l)
        @test L[1] != 1
        @test all(i->i ∈ L, [8, 1, 4, 9])

        cplx, w = witness(X, l, 0.005, false, firstpoint = 1)
        @test size(cplx, 0) == l
        @test size(cplx, 1) == 4
        @test w === nothing

        cplx, w = witness(X, l, 0.1, false, firstpoint = 1)
        @test size(cplx, 0) == l
        @test size(cplx, 1) == 5
        @test w === nothing

        cplx, w = witness(X, l, 0.5, firstpoint = 1)
        @test size(cplx, 0) == l
        @test size(cplx, 1) == 6
        @test size(cplx, 2) == 4
        @test w[0][1] == 0.
        @test w[1][1] ≈ 0.22805075942971
        @test w[2][end-1] ≈ 0.0782128253900049

        cplx2, w2 = witness(X, l, 0.5, expansion=:inductive, firstpoint = 1)
        @test size(cplx2, 0) == l
        @test size(cplx2, 1) == 6
        @test size(cplx2, 2) == 4
        @test w2[0][1] == 0.
        @test w2[1][1] ≈ 0.22805075942971
        @test w2[2][end-1] ≈ 0.0782128253900049

        @test_throws ArgumentError witness(X, l, 0.4, expansion=:a)
    end
end

