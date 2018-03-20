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

    # test iterator
    @testset "Iterator" for (v,c) in flt
        if v == 1
            @test size(c) == (1,)
        elseif  v == 2
            @test size(c) == (2,)
        elseif  v == 3
            @test size(c) == (2,1)
        else
            @test size(c) == (3,2)
        end
    end

    # test io
    let io = IOBuffer()
        write(io, flt)
        @test String(io) == "1,1\n2,2\n1,2,3\n1,3,4\n3,4\n"
        seekstart(io)
        tmp = read(io, Filtration{SimplicialComplex{Int}, Int})
        @testset for ((v1,c1),(v2,c2)) in zip(flt, tmp)
            @test v1 == v2
            @test c1.cells == c2.cells
        end
    end

    # compute boundary matrix
    ∂ = boundary_matrix(flt, reduced=false)
    @test length(∂) == sum(size(flt.complex))
    @test countnz(sparse(∂)) == 4

    # create complex
    cplx = SimplicialComplex(Char)
    push!(cplx, Simplex('a','b','c'), recursive=true)

    # create filtration from an existed complex
    flt = filtration(cplx, Int)

    # compute boundary matrix
    ∂ = boundary_matrix(flt)
    @test length(∂) == 8
    @test countnz(sparse(∂)) == 12

    srand(9236493643764)
    N = 10
    X = rand(3,N)
    cplx, w = vietorisrips(X, 0.4, true)
    flt = filtration(cplx, w)
    ∂ = boundary_matrix(flt)
    @test length(∂) == sum(size(cplx))+1
    @test countnz(sparse(∂)) == 29
end
