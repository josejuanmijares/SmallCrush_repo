include("init_test.jl")
g=RNGTest.wrap(SuperJuice,true,Float64);
println(RNGTest.svaria_WeightDistrib(g, 1, 200000, 27, 256, 0.0, 0.125))

