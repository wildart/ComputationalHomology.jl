@testset "Simplex" begin
    s = Simplex(1, 2, 3)
    @test s == Simplex(1, 2, 3)
    @test sort!(values(s)) == [1,2,3]
    @test eltype(s) == Int
    @test eltype(Simplex{1,Float64}) == Float64

    s = Simplex([1,2,3])
    @test dim(s) == 2
    @test sort!(values(s)) == [1,2,3]
    allfaces = [Simplex(2, 3), Simplex(1, 3), Simplex(1, 2)]
    @testset "Simplex faces" for o in faces(s)
        idx = findfirst(isequal(o), allfaces)
        @test values(o) == values(allfaces[idx])
    end

    ss = Set{Simplex}()
    push!(ss, s)
    push!(ss, s)
    @test length(ss) == 1

    X = hcat(zeros(3), Matrix(I,3,3))
    S = Simplex(1,2,3,4)
    @test volume(S,X) == 1/6

    s = Simplex(1,3)
    fs = collect(faces(s))
    @test dim(s) == 1
    @test all(f -> f ∈ [Simplex(1),Simplex(3)], fs)

    @test fs[2] ∪ fs[1]  == s

end

@testset "Cube" begin
    c = Cube{Vector{Int}}([0,0,0],Int[]) # 0-cube
    @test dim(c) == 0
    c = Cube{Vector{Int}}([0,0,0],Int[0,0,0]) # 0-cube
    @test dim(c) == 0
    c = Cube{Vector{Int}}([0,0,0],[1]) # 1-cube
    @test dim(c) == 1
    c = Cube{Vector{Int}}([0,0,0],[1,0,1]) # 2-cube
    @test dim(c) == 2
    c = Cube{Vector{Int}}([0,0,0],[1,1,1]) # 3-cube
    @test dim(c) == 3
    @test hash(c) == hash([0,0,0,1,1,1])
    @test c == Cube{Vector{Int}}([0,0,0],[1,1,1])

    c = Cube{Vector{Float64}}([0.,0.,0.],[.5,.0,.5]) # 2-cube
    @test volume(c) == 0.5^2
    c = Cube{Vector{Float64}}([0.,0.,0.],[.5,.5,.5]) # 2-cube
    @test volume(c) == 0.5^3

end

@testset "CW" begin
    ComputationalHomology.resetCellId!()
    a = Cell()
    b = Cell()
    @test dim(a) == 0
    @test hash(a) == 0x02011ce34bce797f
    @test a == a
    @test eltype(a) == Int
    @test eltype(Cell{Float64}) == Float64

    @test_throws AssertionError Cell(0, a)

    c = Cell(1, a, b)
    @test dim(c) == 1
    @test vertices(c)[1] == b
    @test vertices(c)[2] == a
    cb = boundary(c)
    @test cb[hash(b)] == -1
    @test cb[hash(a)] == 1

    c2 = Cell()
    d = Cell(1, a, b, c2)
    e = c ∪ d
    @test dim(e) == 2

    @testset for (f1, f2) in zip(vertices(c), [b, a])
        @test f1 == f2
    end

    @testset for (f1, f2) in zip(faces(e), [c, d])
        @test f1 == f2
    end

    @test vertices(e)[1] == c2
    @test vertices(e)[2] == b
    @test vertices(e)[3] == a

    f = Cell(1, [a, b], [1, 1])
    cb = boundary(f)
    @test cb[hash(a)] == 1
    @test cb[hash(b)] == 1
end
