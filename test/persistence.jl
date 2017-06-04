# include("src/ComputationalHomology.jl")
# using ComputationalHomology
# using Base.Test
# include("test\\runtests.jl")
# include("test/runtests.jl")

@testset "Persistent Homology" begin

    # Create complex and filtration simultaneously
    flt = Filtration(SimplicialComplex{Int}, Int)
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
    @test size(flt.complex) == (3,3,1)

    # Compute persistence boundary matrix
    ∂ = boundary_matrix(flt, reduced=false)

    # Twist reduction
    ps, R = pairs(TwistReduction, deepcopy(∂))
    @testset for (p, t) = zip(ps, [2=>4, 3=>5, 6=>7])
        @test p == t
    end

    # Standard reduction
    ps, R = pairs(StandardReduction, deepcopy(∂))
    @testset for (p, t) = zip(ps, [2=>4, 3=>5, 6=>7])
        @test p == t
    end

    flt = Filtration(SimplicialComplex{Int}, Int)
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

    ps, R = pairs(TwistReduction, boundary_matrix(flt, reduced=false))
    itr = ComputationalHomology.endpoints(flt, ps)
    @test itr[0] == [1=>2]
    @test itr[1] == [3=>7]


    flt = Filtration(SimplicialComplex{Int}, Int)
    splxs=[ Simplex(1) => 1,
            Simplex(2) => 2,
            Simplex(3) => 3,
            Simplex(4) => 4,
            Simplex(5) => 5,
            Simplex(1,2) => 6,
            Simplex(2,3) => 7,
            Simplex(3,4) => 8,
            Simplex(4,1) => 9,
            Simplex(3,5) => 10,
            Simplex(4,5) => 11,
            Simplex(3,4,5) => 12]
    for s in splxs
        push!(flt, s...)
    end
    @test size(flt.complex) == (5,6,1)

    ∂ = boundary_matrix(flt)
    R = reduce(StandardReduction, deepcopy(∂))
    @test ComputationalHomology.betti(∂, R, 0) == 1
    @test ComputationalHomology.betti(∂, R, 1) == 1
    @test ComputationalHomology.betti(∂, R, 2) == 0

    ph = persistenthomology(flt, TwistReduction)
    @test group(ph, 0) == 1
    @test group(ph, 1) == 1
    @test group(ph, 2) == 0

    @testset "Method Comparison" for (g1, g2) in zip(homology(flt.complex), persistenthomology(flt, TwistReduction))
        @test g1[1] == g2[1]
        @test g1[2] == g2[2]
    end
end
