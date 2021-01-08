@testset "Tropic Coordinates" begin
    # example from arXiv:1709.02647, p.13
    pd1 = Dict(0 => diagram(1.0=>3.0, 3.0=>4.0))
    pd2 = Dict(0 => diagram(2.0=>4.0))
    pds = map(ComputationalHomology.flatten, [pd1, pd2])
    d,n,m,mv = ComputationalHomology.tropicstats(pds)
    @test m == 3
    @test tropic(pds[1]) == [2, 3, 4, 6, 7]
    xs = tropic(pds)
    @test xs[:,1] == [2, 3, 4, 6, 7]
    @test xs[:,2] == [2, 2, 4, 4, 4]

    # with half-open interval
    img = reshape([0,0,2,2,2,3,0,0,0,2,0,0,1,3,0,0,0,0,2,1,2,0,0,0,0,0,0,4,1,0,0,0,0,0,1,1,1,1,0,0,0,0,2,0,0,2,1,0,0,0,2,0,0,0,2,0,0,0,1,2,2,3,0,0],8,8)
    flt = filtration(BitmapComplex(maximum(img) .- img))
    pdn = diagram(flt)
    @test tropic(pdn) == [4, 6, 4, 8, 8]

    # img = 255 .- readdlm("misc/du8.csv", ',', UInt8)
    # x = BitmapComplex(img) |> filtration |> diagram |> tropic
end
