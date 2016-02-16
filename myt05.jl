include("init_test.jl")
g=RNGTest.wrap(SuperJuice,true,Float64);
println(RNGTest.sknuth_CouponCollector(g, 1, 500000, 26, 16))

