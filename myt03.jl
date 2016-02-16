include("init_test.jl")
g=RNGTest.wrap(SuperJuice,true,Float64);
println(RNGTest.sknuth_Gap(g, 1, 200000, 22, 0.0, .00390625))

