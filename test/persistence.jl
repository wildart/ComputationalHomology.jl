@testset "Persistence" begin

    # Intervals
    p = Interval(1.0 => 2.0)
    q = Interval(0.0 => 3.0)
    r = Interval(0.0 => 4.0)

    @test q < p
    @test q < r

    @test ComputationalHomology.birth(p) == -1.
    @test ComputationalHomology.death(p) == 3.
    @test diag(p)  == Interval(1.5 => 1.5)

    #######
    # from "Computational Topology - An Introduction" by Edelsbrunner & Harer, p. 184
    #######
    flt = Filtration(SimplicialComplex{Simplex{Int}}, Int)
    splxs=[ Simplex(1) => 1,
            Simplex(2) => 2,
            Simplex(3) => 3,
            Simplex(1,2) => 4,
            Simplex(2,3) => 5,
            Simplex(3,1) => 6,
            Simplex(1,2,3) => 7]
    for s in splxs
        push!(flt, s...)
    end
    @test size(complex(flt)) == (3,3,1)

    # Twist reduction
    ps, R = ComputationalHomology.pairs(TwistReduction, flt, reduced=true)
    @testset for (p, t) = zip(ps, [0=>1, 2=>4, 3=>5, 6=>7])
        @test p == t
    end

    # Standard reduction
    ps, R = ComputationalHomology.pairs(StandardReduction, flt)
    @testset for (p, t) = zip(ps, [2=>4, 3=>5, 6=>7, 1=>Inf])
        @test p == t
    end

    #######
    # from "Topology for Computing" by Zomorodian, pp.138-145
    #######
    flt = Filtration(SimplicialComplex{Simplex{Char}}, Int)
    splxs=[ Simplex('a') => 0,
            Simplex('b') => 0,
            Simplex('c') => 1,
            Simplex('d') => 1,
            Simplex('a', 'b') => 1,
            Simplex('b', 'c') => 1,
            Simplex('c', 'd') => 2,
            Simplex('a', 'd') => 2,
            Simplex('a', 'c') => 3,
            Simplex('a', 'b', 'c') => 4,
            Simplex('a', 'c', 'd') => 5]
    for s in splxs
        push!(flt, s...)
    end

    @testset "Intervals " for (d, itrs) in intervals(flt, length0=true)
        titrs = d == 0 ? intervals(0, 0=>Inf, 0=>1, 1=>1, 1=>2) : intervals(1, 3=>4, 2=>5)
        for itr in itrs
            @test itr ∈ titrs
        end
    end

    #######
    # Ex.1
    #######
    flt = Filtration(SimplicialComplex{Simplex{Int}}, Int)
    splxs=[ Simplex(1) => 0,
            Simplex(2) => 0,
            Simplex(3) => 0,
            Simplex(4) => 0,
            Simplex(5) => 1,
            Simplex(1,2) => 0,
            Simplex(2,3) => 0,
            Simplex(3,4) => 0,
            Simplex(4,1) => 0,
            Simplex(3,5) => 2,
            Simplex(4,5) => 3,
            Simplex(3,4,5) => 7]
    for s in splxs
        push!(flt, s...)
    end

    itr = intervals(flt, length0=true)
    @test Interval(0=>0) ∈ itr[0]
    @test Interval(1=>2) ∈ itr[0]
    @test Interval(1, 3=>7) ∈ itr[1]

    itr = intervals(flt, reduction = StandardReduction)
    @test Interval(1=>2) ∈ itr[0]
    @test Interval(1, 3=>7) ∈ itr[1]
    @test Interval(1, 0=>Inf) ∈ itr[1]

    #######
    # Ex.2
    #######
    flt = Filtration(SimplicialComplex{Simplex{Int}}, Float64)
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
    @testset "Method Comparison" for (g1, g2) in zip(homology(complex(flt), Int), persistenthomology(TwistReduction, flt))
        @test g1[1] == g2[1]
        @test g1[2] == g2[2]
    end
end
