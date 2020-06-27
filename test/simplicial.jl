@testset "Simplicial Complex" begin
    cplx = SimplicialComplex()
    @test eltype(cplx) == Simplex
    @test length(cells(cplx)) == 0
    @test length(cells(cplx,2)) == 0
    @test size(cplx) == ()

    cplx = SimplicialComplex(Simplex(1,2,3), Simplex(2,4), Simplex(4))

    @test eltype(cplx) == Simplex
    @test length(cplx) == 9
    @test length(cells(cplx)) == 3
    @test length(cells(cplx,2)) == 1
    @test length(cells(cplx,3)) == 0
    @test size(cplx) == (4,4,1)
    @test size(cplx,-1) == 0
    @test size(cplx,0) == 4
    @test size(cplx,1) == 4
    @test size(cplx,2) == 1
    @test size(cplx,3) == 0
    @test Simplex(4) ∈ cplx
    @test ComputationalHomology.position(cplx, Simplex(4)) == 4
    @test Simplex(10) ∉ cplx

    @test first(Iterators.drop(cplx, 5)) == Simplex(1,3)

    cplx = SimplicialComplex(Simplex(1,2,3),
                       Simplex(2,4),
                       Simplex(3,4),
                       Simplex(5,4),
                       Simplex(6))

    # boundary matrix
    @test convert(Matrix, boundary(cplx, 2, Int)) == [1 -1 1 0 0 0]'

    # simplex values
    cplx = SimplicialComplex()

    splx = push!(cplx, Simplex(1))
    @test ComputationalHomology.position(cplx, splx[1]) == 1
    @test dim(splx[1]) == 0
    @test sum(size(cplx)) == 1
    @test values(cells(cplx, 0)[1]) == [1]

    splx = push!(cplx, Simplex(1,2))
    @test splx[1] == Simplex(1,2)
    @test dim(splx[1]) == 1
    @test sum(size(cplx)) == 2
    @test sort!(values(cells(cplx, 1)[1])) == [1,2]

    cplx = SimplicialComplex()
    splxs = push!(cplx, Simplex(1, 2), recursive=true)
    @test splxs[1] == Simplex(1,2)
    @test sort!(values(splxs[1])) == [1, 2]
    @test splxs[3] == Simplex(2)
    @test values(splxs[3]) == [2]
    @test splxs[2] == Simplex(1)
    @test values(splxs[2]) == [1]
    @test sum(size(cplx)) == 3
    @test size(cplx, 0) == 2
    @test size(cplx, 1) == 1

    splxs = push!(cplx, Simplex(1,3), recursive=true)
    @test splxs[1] == Simplex(1,3)
    @test sort!(values(splxs[1])) == [1,3]
    @test splxs[2] == Simplex(3)
    @test values(splxs[2]) == [3]
    @test sum(size(cplx)) == 5
    @test size(cplx, 0) == 3
    @test size(cplx, 1) == 2

    @test Simplex(1,3) in cplx
    @test first(cells(cplx,1)) in cplx
    @test !(Simplex(2,3) in cplx)
    @test !(Simplex(2,3,3) in cplx)

    io = IOBuffer()
    write(io, cplx)
    seekstart(io)
    cplx2 = read(io, SimplicialComplex(), Simplex{0,Int})
    @test size(cplx2, 0) == 3
    @test size(cplx2, 1) == 2
    @test all(cells(cplx2,1) .== cells(cplx,1))

    A = ComputationalHomology.adjacency_matrix(ComputationalHomology.sphere(1), UInt8)
    @test collect(A) == UInt8[0x00 0x01; 0x01 0x00]
    @test eltype(A) == UInt8
end
