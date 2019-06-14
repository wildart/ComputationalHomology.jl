@testset "CW Complex" begin

    cplx = CWComplex()
    @test size(cplx) == ()
    @test celltype(cplx) == Cell{Int}

    v1 = Cell()
    v2 = Cell()
    cplx = CWComplex(v1, v2)
    @test size(cplx) == (2,)
    @test cells(cplx) == [[v1, v2]]
    @test cells(cplx, 0) == [v1, v2]
    @test celltype(cplx) == Cell{Int}

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

    flt = Filtration(CWComplex, Int)
    e3 = Cell(1, v1, v2)
    e4 = Cell(1, v1, v2)
    f5 = Cell(2, e3, e4)
    f6 = Cell(2, e3, e4)
    push!(flt, v1, 1)
    push!(flt, v2, 2)
    push!(flt, e3, 3)
    push!(flt, e4, 4)
    push!(flt, f5, 5)
    push!(flt, f6, 6)

    @test size(complex(flt)) == (2,2,2)
    homitr = Dict( 0 => intervals(0, 1=>Inf, 2=>3), 1 => intervals(1, 4=>5), 2 => intervals(2, 6=>Inf) )
    @testset "Intervals " for (d, itrs) in intervals(flt)
        for itr in itrs
            @test itr ∈ homitr[d]
        end
    end
end
