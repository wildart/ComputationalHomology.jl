@testset "CW Complex" begin

    cplx = CWComplex()
    @test size(cplx) == ()
    @test eltype(cplx) == Cell

    v1 = Cell()
    v2 = Cell()
    cplx = CWComplex(v1, v2)
    @test size(cplx) == (2,)
    @test cells(cplx) == [[v1, v2]]
    @test cells(cplx, 0) == [v1, v2]
    @test eltype(cplx) == Cell
    @testset "Iterator" for (i,c) in enumerate(cplx)
        @test c == (i == 1 ? v1 : v2)
    end

    cplx = CWComplex()
    push!(cplx, v1)
    push!(cplx, v2)
    @test size(cplx) == (2,)
    @test cells(cplx) == [[v1, v2]]
    @test cells(cplx, 0) == [v1, v2]

    cplx = CWComplex()
    v = push!(cplx, Cell())
    e = push!(cplx, Cell(1, v, [0]))
    @test size(cplx) == (1,1)
    @test cells(cplx, 0) == v
    @test cells(cplx, 1) == e

    flt = Filtration(CWComplex)
    e3 = Cell(1, v1, v2)
    e4 = Cell(1, v1, v2)
    f5 = Cell(2, e3, e4)
    f6 = Cell(2, e3, e4)
    push!(flt, v1, 1.0)
    push!(flt, v2, 2.0)
    push!(flt, e3, 3.0)
    push!(flt, e4, 4.0)
    push!(flt, f5, 5.0)
    push!(flt, f6, 6.0)

    cplx = complex(flt)
    @test size(cplx) == (2,2,2)

    # (co)faces
    @testset "Faces" for fi in faces(cplx, hash(f5))
        @test fi ∈ map(hash, [e3, e4])
    end
    @test length(faces(cplx, hash(v1))) == 0
    @testset "Cofaces" for fi in cofaces(cplx, hash(e4))
        @test fi ∈ map(hash, [f5, f6])
    end
    @test length(cofaces(cplx, hash(f6))) == 0

    homitr = Dict( 0 => diagram(1.0=>Inf, 2.0=>3.0), 1 => diagram(4.0=>5.0), 2 => diagram(6.0=>Inf) )
    @testset "Intervals " for (d, itrs) in diagram(flt)
        for itr in itrs
            @test itr ∈ homitr[d]
        end
    end
end
