#push!(LOAD_PATH, pwd()); 
#LOAD_PATH=unique(LOAD_PATH);
@everywhere push!(LOAD_PATH, pwd()); 
@everywhere LOAD_PATH=unique(LOAD_PATH);
@everywhere using myGA4
@everywhere using RNGTest
@everywhere using Base.Test
