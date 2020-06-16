using ComputationalHomology
using BenchmarkTools
using Random

Random.seed!(2903872398473);
rows, cols = 2, 15
Xsml = rand(rows, cols)
Xmed = rand(rows, cols*10)
Xbig = rand(rows, cols*100)
l = 20
ɛ = 0.05
Csml, wsml = vietorisrips(Xsml, ɛ, maxoutdim=1)
Cmed, wmed = vietorisrips(Xmed, ɛ, maxoutdim=1)
Cbig, wbig = vietorisrips(Xbig, ɛ, maxoutdim=1)

suite = BenchmarkGroup()

suite["parts"] = BenchmarkGroup(["complex"])
suite["parts"]["position-small-start"] = @benchmarkable ComputationalHomology.position(Csml, Csml.cells[0][1])
suite["parts"]["position-small-center"] = @benchmarkable ComputationalHomology.position(Csml, Csml.cells[0][length(Csml.cells[0])>>1])
suite["parts"]["position-small-end"] = @benchmarkable ComputationalHomology.position(Csml, Csml.cells[0][length(Csml.cells[0])])
suite["parts"]["position-medium-start"] = @benchmarkable ComputationalHomology.position(Cmed, Cmed.cells[0][1])
suite["parts"]["position-medium-center"] = @benchmarkable ComputationalHomology.position(Cmed, Cmed.cells[0][length(Cmed.cells[0])>>1])
suite["parts"]["position-medium-end"] = @benchmarkable ComputationalHomology.position(Cmed, Cmed.cells[0][length(Cmed.cells[0])])
suite["parts"]["position-big-start"] = @benchmarkable ComputationalHomology.position(Cbig, Cbig.cells[0][1])
suite["parts"]["position-big-center"] = @benchmarkable ComputationalHomology.position(Cbig, Cbig.cells[0][length(Cbig.cells[0])>>1])
suite["parts"]["position-big-end"] = @benchmarkable ComputationalHomology.position(Cbig, Cbig.cells[0][length(Cbig.cells[0])])

ɛ = 0.25
Asml = ComputationalHomology.adjacency_matrix(Csml, Bool)
Amed = ComputationalHomology.adjacency_matrix(Cmed, Bool)
Abig = ComputationalHomology.adjacency_matrix(Cbig, Bool)

suite["construction"] = BenchmarkGroup(["vrips", "witness"])
suite["construction"]["lowernbrs-small"] = @benchmarkable ComputationalHomology.lowernbrs(Csml, hash(Csml.cells[0][10]), Asml)
suite["construction"]["lowernbrs-medium"] = @benchmarkable ComputationalHomology.lowernbrs(Cmed, hash(Cmed.cells[0][103]), Amed)
suite["construction"]["lowernbrs-big"] = @benchmarkable ComputationalHomology.lowernbrs(Cbig, hash(Cbig.cells[0][1009]), Abig)
suite["construction"]["vr-dim1"] = @benchmarkable vietorisrips(Xmed, ɛ, false, maxoutdim=1)
suite["construction"]["vr-incremental"] = @benchmarkable vietorisrips(Xmed, ɛ, false, maxoutdim=2, expansion=:incremental)
suite["construction"]["vr-inductive"] = @benchmarkable vietorisrips(Xmed, ɛ, false, maxoutdim=2, expansion=:inductive)
suite["construction"]["witness-dim1"] = @benchmarkable witness(Xmed, l, ɛ, false, maxoutdim=1)
suite["construction"]["witness-incremental"] = @benchmarkable witness(Xmed, l, ɛ, false, maxoutdim=2, expansion=:incremental)
suite["construction"]["witness-inductive"] = @benchmarkable witness(Xmed, l, ɛ, false, maxoutdim=2, expansion=:inductive)

Fsml = filtration(Csml, wsml)
Fmed = filtration(Cmed, wmed)
Fbig = filtration(Cbig, wbig)

suite["persistence"] = BenchmarkGroup(["homology", "diagram"])
suite["persistence"]["standard-med"] = @benchmarkable diagram(StandardReduction, Fmed)
suite["persistence"]["standard-big"] = @benchmarkable diagram(StandardReduction, Fbig)
suite["persistence"]["twist-med"] = @benchmarkable diagram(TwistReduction, Fmed)
suite["persistence"]["twist-big"] = @benchmarkable diagram(TwistReduction, Fbig)


# parameters
loadparams!(suite, BenchmarkTools.load("bmsetup.json")[1], :evals, :samples);
# tune!(suite);
# BenchmarkTools.save("bmsetup.json", params(suite));

# benchmark
results = run(suite, verbose = true)
BenchmarkTools.save("results.json", results)
