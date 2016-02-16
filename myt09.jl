include("init_test.jl")
g=RNGTest.wrap(SuperJuice,true,Float64);
println(RNGTest.sstring_HammingIndep(g, 1, 500000, 20, 10, 300, 0)[:Mean])

