include("init_test.jl")
g=RNGTest.wrap(SuperJuice,true,Float64);
println(RNGTest.sknuth_SimpPoker(g, 1, 400000, 24, 64, 64))

