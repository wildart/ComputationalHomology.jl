@testset "Filtration" begin
    # create empty filtration
    flt = Filtration(SimplicialComplex{Int}, Int)

    # fill it with (cell, 'filtration value') pairs
    push!(flt, Simplex(1), 1)
    @test size(complex(flt)) == (1,)
    push!(flt, Simplex(2), 2)
    @test size(complex(flt)) == (2,)
    push!(flt, Simplex(1,2), 3, recursive=true)
    @test size(complex(flt)) == (2,1)
    push!(flt, Simplex(1,3), 4, recursive=true)
    @test size(complex(flt)) == (3,2)

    # test iterator
    @testset "Iterator" for (v,c) in flt
        if v == 1
            @test size(c) == (1,)
        elseif  v == 2
            @test size(c) == (2,)
        elseif  v == 3
            @test size(c) == (2,1)
        else
            @test sum(size(c)) > 3
        end
    end

    # test io
    let io = IOBuffer()
        write(io, flt)
        @test String(take!(copy(io))) == "1,1\n2,2\n1,2,3\n3,4\n1,3,4\n"
        seekstart(io)
        tmp = read(io, Filtration{SimplicialComplex{Int}, Int})
        @testset for ((v1,c1),(v2,c2)) in zip(flt, tmp)
            @test v1 == v2
            @test c1.cells == c2.cells
        end
    end

    # compute boundary matrix
    ∂ = boundary_matrix(flt)
    @test length(∂) == sum(size(complex(flt)))
    @test count(!iszero, sparse(∂)) == 4

    # create filtration from an existed complex
    cplx = SimplicialComplex(Char)
    push!(cplx, Simplex('a','b','c'), recursive=true)
    flt = filtration(cplx)
    @test flt.total[end] == (2, 1, 7)

    # compute boundary matrix
    ∂ = boundary_matrix(flt)
    @test length(∂) == 7
    @test count(!iszero, sparse(∂)) == 9

    Random.seed!(9236493643764)
    N = 10
    X = rand(3,N)
    cplx, w = vietorisrips(X, 0.4, true)
    flt = filtration(cplx, w)
    ∂ = boundary_matrix(flt, reduced=true)
    @test length(∂) == sum(size(cplx))+1
    @test count(!iszero, sparse(∂)) == 29

    SF = simplices(flt)
    @test length(SF) == 9
    @test length(flt) == 19
    @test length((first(v) for v in flt)) == 19
    @test sum(map(e->length(e[2]), SF)) == 19

end
