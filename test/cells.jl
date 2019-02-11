@testset "Simplex" begin
    s = Simplex(1, 2, 3)
    @test s == Simplex(1, 2, 3)
    @test sort!(values(s)) == [1,2,3]

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

    s = Simplex([1,2], [1,3])
    @test dim(s) == 1

    X = hcat(zeros(3), Matrix(I,3,3))
    S = Simplex(1,2,3,4)
    @test volume(S,X) == 1/6

    @test faces(s)[1] ∪ faces(s)[2] == s
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
    @test values(c) == [[0,0,0],[1,1,1]]

    c = Cube{Vector{Float64}}([0.,0.,0.],[.5,.0,.5]) # 2-cube
    @test volume(c) == 0.5^2
    c = Cube{Vector{Float64}}([0.,0.,0.],[.5,.5,.5]) # 2-cube
    @test volume(c) == 0.5^3
end
