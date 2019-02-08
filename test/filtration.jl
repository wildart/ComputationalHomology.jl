@testset "Filtration" begin
    # create empty filtration
    flt = Filtration(SimplicialComplex{Simplex{Int}}, Int)

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
    dimhash(s) = (dim(s), hash(s))
    @testset "Iterator" for ((f1,ss1),(f2,ss2)) in zip(flt,
            enumerate([
                [dimhash(Simplex(1))],
                [dimhash(Simplex(2))],
                [dimhash(Simplex(1,2))],
                [dimhash(Simplex(3)), dimhash(Simplex(1,3))]
            ]))
        @test f1 == f2
        for (s1, s2) in zip(ss1, ss2)
            @test s1 == s2
        end
    end

    # test io
    let io = IOBuffer()
        write(io, flt)
        @test String(take!(copy(io))) == "1,1\n2,2\n2,1,3\n3,4\n3,1,4\n"
        seekstart(io)
        tmp = read(io, Filtration{SimplicialComplex{Simplex{Int}}, Int})
        @testset for ((f1,ss1),(f2,ss2)) in zip(flt, tmp)
            @test f1 == f2
            for (s1, s2) in zip(ss1, ss2)
                @test s1 == s2
            end
        end
    end

    # compute boundary matrix
    ∂ = boundary(flt)
    @test length(∂) == sum(size(complex(flt)))
    @test count(!iszero, sparse(∂)) == 4

    # create filtration from an existed complex
    cplx = SimplicialComplex(Simplex{Char})
    push!(cplx, Simplex('a','b','c'), recursive=true)
    flt = filtration(cplx)
    @test order(flt)[end] == (2, hash(Simplex('a','b','c')), 7)

    # compute boundary matrix
    ∂ = boundary(flt)
    @test length(∂) == 7
    @test count(!iszero, sparse(∂)) == 9

    Random.seed!(9236493643764)
    N = 10
    X = rand(3,N)
    cplx, w = vietorisrips(X, 0.4, true)
    flt = filtration(cplx, w)
    ∂ = boundary(flt)
    @test length(∂) == sum(size(cplx))
    @test count(!iszero, sparse(∂)) == sum(size(cplx))

    SF = simplices(flt)
    @test length(SF) == 9
    @test length(flt) == 9
    @test length(complex(flt)) == 19
    let (fv, ss) = first(flt)
        @test fv == 0.0
        @test length(ss) == N
    end
    @test sum(length(last(v)) for v in flt) == 19
    @test sum(map(e->length(e[2]), SF)) == 19
    @testset for (l1, l2) in zip((length(last(v)) for v in flt), [10, 1, 1, 1, 1, 1, 1, 2, 1])
        @test l1 == l2
    end

    @test isinf(flt.divisions)
    @test minimum(flt) == 0.0
    @test maximum(flt) == mapreduce(e->e[3], max, order(flt))
    flt.divisions = 3
    @test sum(length(last(v)) for v in flt) == 19
    @testset for (l1, l2) in zip((length(last(v)) for v in flt), [11, 2, 6])
        @test l1 == l2
    end
end
