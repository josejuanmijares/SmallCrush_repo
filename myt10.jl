include("init_test.jl")
g=RNGTest.wrap(SuperJuice,true,Float64);
println(RNGTest.swalk_RandomWalk1(g, 1, 1000000, 0, 30, 150, 150))

