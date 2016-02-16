include("init_test.jl")
g=RNGTest.wrap(SuperJuice,true,Float64);
println(RNGTest.sknuth_MaxOft(g, 1, 2000000, 0, 100000, 6))

