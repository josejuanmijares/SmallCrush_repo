include("init_test.jl")

g=RNGTest.wrap(SuperJuice,true,Float64);
println(RNGTest.smarsa_BirthdaySpacings(g, 1, 5000000, 0, 1073741824, 2, 1))
