@testset "Persistence" begin

    # Intervals
    p = Interval(1.0 => 2.0)
    q = Interval(0.0 => 3.0)
    r = Interval(0.0 => 4.0)

    @test q < p
    @test q < r

    @test ComputationalHomology.birthx(p) == -1.
    @test ComputationalHomology.deathx(p) == 3.
    @test diag(p)  == Interval(1.5 => 1.5)

    #######
    # from "Computational Topology - An Introduction" by Edelsbrunner & Harer, p. 184
    #######
    flt = Filtration(SimplicialComplex, Float64)
    splxs=[ Simplex(1) => 1.0,
            Simplex(2) => 2.0,
            Simplex(3) => 3.0,
            Simplex(1,2) => 4.0,
            Simplex(2,3) => 5.0,
            Simplex(3,1) => 6.0,
            Simplex(1,2,3) => 7.0]
    for s in splxs
        push!(flt, s...)
    end
    @test size(complex(flt)) == (3,3,1)

    dgm = Dict( 0 => diagram(2.0=>4.0, 3.0=>5.0, 1.0=>Inf), 1 => diagram(6.0=>7.0))

    @testset "Intervals for TwistReduction" for (d, itrs) in diagram(TwistReduction, flt)
        for itr in itrs
            @test itr ∈ dgm[d]
        end
    end

    @testset "Intervals for Standard Reduction" for (d, itrs) in diagram(StandardReduction, flt)
        for itr in itrs
            @test itr ∈ dgm[d]
        end
    end

    @testset "Intervals for Standard Reduction" for (d, itrs) in diagram(PersistentCocycleReduction{Float64}, flt)
        for itr in itrs
            @test itr ∈ dgm[d]
        end
    end

    #######
    # from "Topology for Computing" by Zomorodian, pp.138-145
    #######
    flt = Filtration(SimplicialComplex)
    splxs=[ Simplex('a') => 0.0,
            Simplex('b') => 0.0,
            Simplex('c') => 1.0,
            Simplex('d') => 1.0,
            Simplex('a', 'b') => 1.0,
            Simplex('b', 'c') => 1.0,
            Simplex('c', 'd') => 2.0,
            Simplex('a', 'd') => 2.0,
            Simplex('a', 'c') => 3.0,
            Simplex('a', 'b', 'c') => 4.0,
            Simplex('a', 'c', 'd') => 5.0]
    for s in splxs
        push!(flt, s...)
    end

    @testset "PH Intervals " for (d, itrs) in diagram(flt, length0=true)
        titrs = d == 0 ? diagram(0.0=>Inf, 0.0=>1.0, 1.0=>1.0, 1.0=>2.0) : diagram(3.0=>4.0, 2.0=>5.0)
        for itr in itrs
            @test itr ∈ titrs
        end
    end

    # PCohomology
    gens = Dict(
        0 => [
                Chain(0, Float64) + (Simplex('a'), 1.0) + (Simplex('b'), 1.0) + (Simplex('c'), 1.0) + (Simplex('d'), 1.0),
                Chain(0, Float64) + (Simplex('b'), 1.0),
                Chain(0, Float64) + (Simplex('d'), 1.0),
            ],
        1 => [
                Chain(1, Float64) + (Simplex('a', 'd'), 1.0),
                Chain(1, Float64) + (Simplex('a', 'c'), 1.0),
            ]
    )
    @testset "PCH Intervals " for ((d1, itrs1), (d2, itrs2)) in zip( diagram(flt), diagram(PersistentCocycleReduction{Float64}, flt))
        for (i1,i2) in zip(
                sort(itrs1, by=ComputationalHomology.birthx),
                sort(itrs2, by=ComputationalHomology.birthx),
            )
            @test i1 == i2
        end
        for (ch1,ch2) in zip(itrs2, gens[d2])
            for (c1, c2) in zip(ch1.g,ch2)
                @test c1 == c2
            end
        end
    end

    #######
    # Ex.1
    #######
    flt = Filtration(SimplicialComplex)
    splxs=[ Simplex(1) => 0.0,
            Simplex(2) => 0.0,
            Simplex(3) => 0.0,
            Simplex(4) => 0.0,
            Simplex(5) => 1.0,
            Simplex(1,2) => 0.0,
            Simplex(2,3) => 0.0,
            Simplex(3,4) => 0.0,
            Simplex(4,1) => 0.0,
            Simplex(3,5) => 2.0,
            Simplex(4,5) => 3.0,
            Simplex(3,4,5) => 7.0]
    for s in splxs
        push!(flt, s...)
    end

    itr = diagram(flt, length0=true)
    @test Interval(0.0=>0.0) ∈ itr[0]
    @test Interval(1.0=>2.0) ∈ itr[0]
    @test Interval(3.0=>7.0) ∈ itr[1]

    itr = diagram(StandardReduction, flt)
    @test Interval(1.0=>2.0) ∈ itr[0]
    @test Interval(3.0=>7.0) ∈ itr[1]
    @test Interval(0.0=>Inf) ∈ itr[1]

    #######
    # Ex.2
    #######
    flt = Filtration(SimplicialComplex)
    splxs=[ Simplex(1) => 1.0,
            Simplex(2) => 2.0,
            Simplex(3) => 3.0,
            Simplex(4) => 4.0,
            Simplex(5) => 5.0,
            Simplex(1,2) => 6.0,
            Simplex(2,3) => 7.0,
            Simplex(3,4) => 8.0,
            Simplex(4,1) => 9.0,
            Simplex(3,5) => 10.0,
            Simplex(4,5) => 11.0,
            Simplex(3,4,5) => 12.0]
    for s in splxs
        push!(flt, s...)
    end
    @test size(complex(flt)) == (5,6,1)
    @test length(flt) == 12

    # Betti numbers
    ∂ = boundary(flt)
    @test sparse(∂)[5,11] == 5
    @test length(∂[9]) > 0
    R = reduce(StandardReduction, ∂)
    @test ComputationalHomology.betti(flt, R, 0) == 1
    @test ComputationalHomology.betti(flt, R, 1) == 1
    @test ComputationalHomology.betti(flt, R, 2) == 0
    @test_throws AssertionError ComputationalHomology.betti(flt, R, 3)
    @test length(R[9]) == 0
    reduce!(StandardReduction, ∂)
    @test length(∂[9]) == 0

    # PH object
    ph = persistenthomology(flt)
    @test eltype(ph) == Tuple{Int64,Int64}
    @test length(ph) == 3
    @test group(ph, 0) == 1
    @test group(ph, 1) == 1
    @test group(ph, 2) == 0
    @test_throws AssertionError group(ph, 3)

    # PH iterator
    @testset "Method Comparison" for (g1, g2) in zip(homology(Int, complex(flt)), persistenthomology(TwistReduction, flt))
        @test g1[1] == g2[1]
        @test g1[2] == g2[2]
    end
end
