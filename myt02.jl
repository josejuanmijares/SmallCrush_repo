include("init_test.jl")
g=RNGTest.wrap(SuperJuice,true,Float64);
println(RNGTest.sknuth_Collision(g, 1, 5000000, 0, 65536, 2))

