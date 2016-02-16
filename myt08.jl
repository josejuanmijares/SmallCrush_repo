include("init_test.jl")
g=RNGTest.wrap(SuperJuice,true,Float64);
println(RNGTest.smarsa_MatrixRank(g, 1, 20000, 20, 10, 60, 60))

