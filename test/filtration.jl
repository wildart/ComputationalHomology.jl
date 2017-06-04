@testset "Filtration" begin

    # create empty filtration
    flt = Filtration(SimplicialComplex{Int}, Int)

    # fill it with (cell, 'filtration value') pairs
    push!(flt, Simplex(1), 1)
    @test size(flt.complex) == (1,)
    push!(flt, Simplex(2), 2)
    @test size(flt.complex) == (2,)
    push!(flt, Simplex(1,2), 3, recursive=true)
    @test size(flt.complex) == (2,1)
    push!(flt, Simplex(1,3), 4, recursive=true)
    @test size(flt.complex) == (3,2)

    # compute boundary matrix
    ∂ = boundary_matrix(flt, reduced=false)
    @test length(∂) == sum(size(flt.complex))
    @test countnz(sparse(∂)) == 4

    # create complex
    cplx = SimplicialComplex(Char)
    push!(cplx, Simplex('a','b','c'), recursive=true)

    # create filtration from an existed complex
    flt = Filtration(cplx, Int)

    # compute boundary matrix
    ∂ = boundary_matrix(flt)
    @test length(∂) == 8
    @test countnz(sparse(∂)) == 12

end
