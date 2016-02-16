module RNGTest

		using Compat
		using myGA4
		
		import Base: convert, getindex, pointer

		const libtestu01 = joinpath(pwd(), "deps", "libtestu01wrapper")
		# const libtestu01 = joinpath(Pkg.dir("RNGTest"), "deps", "libtestu01wrapper")

		swrite = cglobal(("swrite_Basic", libtestu01), Ptr{Bool})
		unsafe_store!(swrite, 0, 1)

		# WrappedRNG

		# TestU01 expects a standard function as input (created via
		# cfunction here). When one wants to test the random stream of
		# an AbstractRNG, a little more work has to be done before
		# passing it to TestU01. The type WrappedRNG wraps an
		# AbstractRNG into an object which knows which type of random
		# numbers to produce, and also whether the scalar or array API
		# should be used (this is useful when different algorithms are
		# used for each case, which results in different streams which
		# should be tested separately). Such an object can then be
		# passed to the Unif01 constructor.

		@compat typealias TestableNumbers Union{Int8, UInt8, Int16, UInt16, Int32, UInt32,
				Int64, UInt64, Int128, UInt128, Float16, Float32, Float64}

		type WrappedRNG{T<:TestableNumbers, RNG<:Any}
				rng::RNG
				cache::Vector{T}
				fillarray::Bool
				vals::Vector{UInt32}
				idx::Int
				TRNGflag::Bool
		end

		function WrappedRNG{RNG, T}(rng::RNG, TRNGflag::Bool, ::Type{T}, fillarray = true, cache_size = 3*2^11 ÷ sizeof(T))
				
				if T <: Integer && cache_size*sizeof(T) % sizeof(UInt32) != 0
						error("cache_size must be a multiple of $(Int(4/sizeof(T))) (for type $T)")
				elseif T === Float16 && cache_size % 6 != 0 || T === Float32 && cache_size % 3 != 0
						error("cache_size must be a multiple of 3 (resp. 6) for Float32 (resp. Float16)")
				end
				cache = Array(T, cache_size)
				
		
				fillcache(WrappedRNG{T, RNG}(rng, cache, fillarray,
							pointer_to_array(convert(Ptr{UInt32}, pointer(cache)), sizeof(cache)÷sizeof(UInt32)),
							0, TRNGflag)) # 0 is a dummy value, which will be set correctly by fillcache
		end
		
		function WrappedRNG{RNG,T}(rng::RNG,rng_f::Function, TRNGflag::Bool, ::Type{T}, fillarray = true, cache_size = 2^14 ÷ sizeof(T))
				
				if T <: Integer && cache_size*sizeof(T) % sizeof(UInt32) != 0
						error("cache_size must be a multiple of $(Int(4/sizeof(T))) (for type $T)")
				elseif T === Float16 && cache_size % 6 != 0 || T === Float32 && cache_size % 3 != 0
						error("cache_size must be a multiple of 3 (resp. 6) for Float32 (resp. Float16)")
				end
				cache = Array(T, cache_size)
				
				fillcache(WrappedRNG{T, RNG}(rng, rng_f, cache, fillarray,
 							pointer_to_array(convert(Ptr{UInt32}, pointer(cache)), sizeof(cache)÷sizeof(UInt32)),
 							0, TRNGflag)) # 0 is a dummy value, which will be set correctly by fillcache
		end
		
		
		

		# The ability to play with the cache size and the fillarray option is for advanced uses,
		# when one wants to test different code path of the particular RNG implementations, like
		# MersenneTwister from Base.
		# For now let's document only the type parameter in the wrap function:
		wrap{T<:TestableNumbers}(rng::AbstractRNG, TRNGflag::Bool, ::Type{T}) = WrappedRNG(rng,TRNGflag, T)
		wrap{T<:TestableNumbers}(rng::Type{myGA4.SuperJuice}, TRNGflag::Bool,::Type{T}) = WrappedRNG(rng,TRNGflag, T)

		function fillcache{T}(g::WrappedRNG{T})
			if g.TRNGflag
				if g.fillarray==true
					#@printf "." 
					g.cache = copy(SuperJuice())
					#gc()
				end
			else
				println("g.TRNGflag is false")
				if g.fillarray
						rand!(g.rng, g.cache)
				else
						for i = 1:length(g.cache)
								@inbounds g.cache[i] = rand(g.rng, T)
						end
				end
			end
			g.idx = 0
			return g
		end

		function call{T<:Integer}(g::WrappedRNG{T})
				g.idx+1 > length(g.vals) && fillcache(g)
				@inbounds return g.vals[g.idx+=1]
		end

		function call(g::WrappedRNG{Float64})
			
				g.idx+1 > length(g.cache) && fillcache(g)
				@inbounds return g.cache[g.idx+=1]
		end

		function call(g::WrappedRNG{Float32})
				g.idx+3 > length(g.cache) && fillcache(g)
				@inbounds begin
						f = Float64(g.cache[g.idx+1])
						# a Float32 has 24 bits of precision, but only 23 bit of entropy
						f += Float64(g.cache[g.idx+2])/exp2(23)
						f += Float64(g.cache[g.idx+=3])/exp2(46)
						return f % 1.0
				end
		end

		function call(g::WrappedRNG{Float16})
				g.idx+6 > length(g.cache) && fillcache(g)
				@inbounds begin
						f = Float64(g.cache[g.idx+1])
						# a Float16 has 10 bits of entropy
						f += Float64(g.cache[g.idx+2])/exp2(10)
						f += Float64(g.cache[g.idx+3])/exp2(20)
						f += Float64(g.cache[g.idx+4])/exp2(30)
						f += Float64(g.cache[g.idx+5])/exp2(40)
						f += Float64(g.cache[g.idx+=6])/exp2(50)
						return f % 1.0
				end
		end


		# Generator type
		type Unif01
				ptr::Ptr{Array{Int32}}
				gentype::Type
				name::ASCIIString
				function Unif01(f::Function, genname)
						for i in 1:100
								tmp = f()
								if typeof(tmp) != Float64 error("Function must return Float64") end
								if tmp < 0 || tmp > 1 error("Function must return values on [0,1]") end
						end
						cf = cfunction(f, Float64, ())
						@compat b = new(ccall((:unif01_CreateExternGen01, libtestu01), Ptr{Void}, (Ptr{UInt8}, Ptr{Void}), genname, cf), Float64)
						#finalizer(b, delete) # TestU01 crashed if two unif01 object are generated. The only safe thing is to explicitly delete the object when used.
						return b
				end

				@compat function Unif01{T<:AbstractFloat}(g::WrappedRNG{T}, genname)
						# we assume that g being created out of an AbstractRNG, it produces Floats in the interval [0,1)
						#@printf "*" 
						@eval f() = call($g) :: Float64
						#@printf "#"
						cf = cfunction(f, Float64, ())
						#@printf "&"
						@compat return new(ccall((:unif01_CreateExternGen01, libtestu01), Ptr{Void}, (Ptr{UInt8}, Ptr{Void}), genname, cf), Float64)
				end
				
				function Unif01{T<:Integer}(g::WrappedRNG{T}, genname)
						@assert Cuint === UInt32
						@eval f() = call($g) :: UInt32
						cf = cfunction(f, UInt32, ())
						@compat return new(ccall((:unif01_CreateExternGenBits, libtestu01), Ptr{Void}, (Ptr{UInt8}, Ptr{Void}), genname, cf), UInt32)
				end
		end
		function delete(obj::Unif01)
						if obj.gentype === Float64
								ccall((:unif01_DeleteExternGen01, libtestu01), Void, (Ptr{Void},), obj.ptr)
						else
								ccall((:unif01_DeleteExternGenBits, libtestu01), Void, (Ptr{Void},), obj.ptr)
						end
		end

		@compat typealias Generator Union{Function, WrappedRNG}

		# Result types

		## gofw
		immutable Gotw_TestArray
				data::Vector{Float64}
				Gotw_TestArray() = new(Array(Float64, 11))
		end
		function getindex(obj::Gotw_TestArray, i::Symbol)
				i == :KSP && return obj.data[1]
				i == :KSM && return obj.data[2]
				i == :KS && return obj.data[3]
				i == :AD && return obj.data[4]
				i == :CM && return obj.data[5]
				i == :WG && return obj.data[6]
				i == :WU	&& return obj.data[7]
				i == :Mean && return obj.data[8]
				i == :Var && return obj.data[9]
				i == :Cor && return obj.data[10]
				i == :Sum && return obj.data[11]
				throw(BoundsError())
		end

		for (t, sCreate, sDelete, sPval) in
				((:ResPoisson, :sres_CreatePoisson, :sres_DeletePoisson, :getPValPoisson),
				 (:MarsaRes, :smarsa_CreateRes, :smarsa_DeleteRes, :getPValSmarsa),
				 (:KnuthRes2, :sknuth_CreateRes2, :sknuth_DeleteRes2, :getPValRes2))
				@eval begin
						# The types
						type $t
								ptr::Ptr{Void}
								function $(t)()
										res = new(ccall(($(string(sCreate)), libtestu01), Ptr{Void}, (), ))
										finalizer(res, delete)
										return res
								end
						end
						# Finalizers
						function delete(obj::$t)
								ccall(($(string(sDelete)), libtestu01), Void, (Ptr{Void},), obj.ptr)
						end
						# pvalue extractors
						pvalue(obj::$t) = ccall(($(string(sPval)), libtestu01), Float64, (Ptr{Void},), obj.ptr)
				end
		end

		# Sres
		# Basic
		## Type
		type ResBasic
				ptr::Ptr{Void}
				function ResBasic()
						res = new(ccall((:sres_CreateBasic, libtestu01), Ptr{Void}, (), ))
						finalizer(res, delete)
						return res
				end
		end
		## Finalizers
		function delete(obj::ResBasic)
				ccall((:sres_DeleteBasic, libtestu01), Void, (Ptr{Void},), obj.ptr)
		end
		## pvalue extractors
		function pvalue(obj::ResBasic)
				res = Gotw_TestArray()
				ccall((:getPValBasic, libtestu01), Void, (Ptr{Void}, Ptr{Float64}), obj.ptr, res.data)
				return res
		end

		# Chi2
		## Type
		type ResChi2
				ptr::Ptr{Void}
				N::Int
				function ResChi2(N::Integer)
						res = new(ccall((:sres_CreateChi2, libtestu01), Ptr{Void}, (), ), N)
						finalizer(res, delete)
						return res
				end
		end
		## Finalizers
		function delete(obj::ResChi2)
				ccall((:sres_DeleteChi2, libtestu01), Void, (Ptr{Void},), obj.ptr)
		end
		## pvalue extractors
		function pvalue(obj::ResChi2)
				res = Gotw_TestArray()
				ccall((:getPValChi2, libtestu01), Void, (Ptr{Void}, Ptr{Float64}), obj.ptr, res.data)
				return res[obj.N == 1 ? :Mean : :Sum]
		end

		# sknuth
		## Type
		type KnuthRes1
				ptr::Ptr{Void}
				N::Int
				function KnuthRes1(N::Integer)
						res = new(ccall((:sknuth_CreateRes1, libtestu01), Ptr{Void}, (), ), N)
						finalizer(res, delete)
						return res
				end
		end
		## Finalizers
		function delete(obj::KnuthRes1)
				ccall((:sknuth_DeleteRes1, libtestu01), Void, (Ptr{Void},), obj.ptr)
		end
		## pvalue extractors
		function pvalue(obj::KnuthRes1)
				chi = Gotw_TestArray()
				bas = Gotw_TestArray()
				ccall((:getPValRes1, libtestu01), Void, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}), obj.ptr, chi.data, bas.data)
				return chi[obj.N == 1 ? :Mean : :Sum], bas[obj.N == 1 ? :Mean : :AD]
		end

		# smarsa
		## Type
		type MarsaRes2
				ptr::Ptr{Void}
				N::Int
				function MarsaRes2(N::Integer)
						res = new(ccall((:smarsa_CreateRes2, libtestu01), Ptr{Void}, (), ), N)
						finalizer(res, delete)
						return res
				end
		end
		## Finalizers
		function delete(obj::MarsaRes2)
				ccall((:smarsa_DeleteRes2, libtestu01), Void, (Ptr{Void},), obj.ptr)
		end
		## pvalue extractors
		function pvalue(obj::MarsaRes2)
				res = Gotw_TestArray()
				ccall((:getPValSmarsa2, libtestu01), Void, (Ptr{Void}, Ptr{Float64}), obj.ptr, res.data)
				return obj.N == 1 ? res[:Mean] : res[:Sum]
		end

		# Walk
		## Type
		type WalkRes
				ptr::Ptr{Void}
				N::Int
				function WalkRes(N::Integer)
						res = new(ccall((:swalk_CreateRes, libtestu01), Ptr{Void}, (), ), N)
						finalizer(res, delete)
						return res
				end
		end
		## Finalizers
		function delete(obj::WalkRes)
				ccall((:swalk_DeleteRes, libtestu01), Void, (Ptr{Void},), obj.ptr)
		end
		## pvalue extractors
		function pvalue(obj::WalkRes)
				pH = Gotw_TestArray()
				pM = Gotw_TestArray()
				pJ = Gotw_TestArray()
				pR = Gotw_TestArray()
				pC = Gotw_TestArray()
				ccall((:getPVal_Walk, libtestu01), Void, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), obj.ptr, pH.data, pM.data, pJ.data, pR.data, pC.data)
				return obj.N == 1 ? (pH[:Mean], pM[:Mean], pJ[:Mean], pR[:Mean], pC[:Mean]) : (pH[:Sum], pM[:Sum], pJ[:Sum], pR[:Sum], pC[:Sum])
		end

		# Npairs
		immutable Snpair_StatArray
				data::Vector{Float64}
				Snpair_StatArray() = new(Array(Float64, 11))
		end
		function getindex(obj::Snpair_StatArray, i::Symbol)
				i == :NP && return obj.data[1]
				i == :NPS && return obj.data[2]
				i == :NPPR && return obj.data[3]
				i == :mNP && return obj.data[4]
				i == :mNP1 && return obj.data[5]
				i == :mNP1S && return obj.data[6]
				i == :mNP2	&& return obj.data[7]
				i == :mNP2S && return obj.data[8]
				i == :NJumps && return obj.data[9]
				i == :BB && return obj.data[10]
				i == :BM && return obj.data[11]
				throw(BoundsError())
		end
		## Type
		type NpairRes
				ptr::Ptr{Void}
				N::Int
				function NpairRes(N::Integer)
						res = new(ccall((:snpair_CreateRes, libtestu01), Ptr{Void}, (), ), N)
						finalizer(res, delete)
						return res
				end
		end
		## Finalizers
		function delete(obj::NpairRes)
				ccall((:snpair_DeleteRes, libtestu01), Void, (Ptr{Void},), obj.ptr)
		end
		## pvalue extractors
		function pvalue(obj::NpairRes)
				res = Snpair_StatArray()
				ccall((:getPVal_Npairs, libtestu01), Void, (Ptr{Void}, Ptr{Float64}), obj.ptr, res.data)
				return obj.N == 1 ? (res[:NP], res[:mNP]) : (res[:NP], res[:mNP1], res[:mNP2], res[:NJumps])
		end

		# scomp
		## Type
		type CompRes
				ptr::Ptr{Void}
				N::Int
				function CompRes(N::Integer)
						res = new(ccall((:scomp_CreateRes, libtestu01), Ptr{Void}, (), ), N)
						finalizer(res, delete)
						return res
				end
		end
		## Finalizers
		function delete(obj::CompRes)
				ccall((:scomp_DeleteRes, libtestu01), Void, (Ptr{Void},), obj.ptr)
		end
		## pvalue extractors
		function pvalue(obj::CompRes)
				num = Gotw_TestArray()
				size = Gotw_TestArray()
				ccall((:getPValScomp, libtestu01), Void, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}), obj.ptr, num.data, size.data)
				return num[obj.N == 1 ? :Mean : :Sum], size[obj.N == 1 ? :Mean : :Sum]
		end

		# sspectral
		## Type
		type SpectralRes
				ptr::Ptr{Void}
				function SpectralRes()
						res = new(ccall((:sspectral_CreateRes, libtestu01), Ptr{Void}, (), ))
						finalizer(res, delete)
						return res
				end
		end
		## Finalizers
		function delete(obj::SpectralRes)
				ccall((:sspectral_DeleteRes, libtestu01), Void, (Ptr{Void},), obj.ptr)
		end
		## pvalue extractors
		function pvalue(obj::SpectralRes)
				res = Gotw_TestArray()
				ccall((:getPValSspectral, libtestu01), Void, (Ptr{Void}, Ptr{Float64}), obj.ptr, res.data)
				return res[:AD]
		end

		# sstring
		## Type
		type StringRes
				ptr::Ptr{Void}
				function StringRes()
						res = new(ccall((:sstring_CreateRes, libtestu01), Ptr{Void}, (), ))
						finalizer(res, delete)
						return res
				end
		end
		## Finalizers
		function delete(obj::StringRes)
				ccall((:sstring_DeleteRes, libtestu01), Void, (Ptr{Void},), obj.ptr)
		end
		## pvalue extractors
		function pvalue(obj::StringRes)
				res = Gotw_TestArray()
				ccall((:getPValStringRes, libtestu01), Void, (Ptr{Void}, Ptr{Float64}), obj.ptr, res.data)
				return res
		end

		## Type
		type StringRes2
				ptr::Ptr{Void}
				N::Int
				function StringRes2(N::Integer)
						res = new(ccall((:sstring_CreateRes2, libtestu01), Ptr{Void}, (), ), N)
						finalizer(res, delete)
						return res
				end
		end
		## Finalizers
		function delete(obj::StringRes2)
				ccall((:sstring_DeleteRes2, libtestu01), Void, (Ptr{Void},), obj.ptr)
		end
		## pvalue extractors
		function pvalue(obj::StringRes2)
				res = Gotw_TestArray()
				res2 = ccall((:getPValStringRes2, libtestu01), Float64, (Ptr{Void}, Ptr{Float64}), obj.ptr, res.data)
				return res[obj.N == 1 ? :Mean : :Sum], res2
		end

		## Type
		type StringRes3
				ptr::Ptr{Void}
				N::Int
				function StringRes3(N::Integer)
						res = new(ccall((:sstring_CreateRes3, libtestu01), Ptr{Void}, (), ), N)
						finalizer(res, delete)
						return res
				end
		end
		## Finalizers
		function delete(obj::StringRes3)
				ccall((:sstring_DeleteRes3, libtestu01), Void, (Ptr{Void},), obj.ptr)
		end
		## pvalue extractors
		function pvalue(obj::StringRes3)
				res1 = Gotw_TestArray()
				res2 = Gotw_TestArray()
				ccall((:getPValStringRes3, libtestu01), Float64, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}), obj.ptr, res1.data, res2.data)
				return res1[obj.N == 1 ? :Mean : :Sum], res2[obj.N == 1 ? :Mean : :Sum]
		end

		#########
		# Tests #
		#########
		## smarsa
		function smarsa_BirthdaySpacings(gen::Generator, N::Integer, n::Integer, r::Integer, d::Integer, t::Integer, p::Integer)
				 @printf "\nTEST smarsa_BirthdaySpacings \t"
				 unif01 = Unif01(gen, "")
				 sres = ResPoisson()
				 ccall((:smarsa_BirthdaySpacings, libtestu01), Void,
						(Ptr{Void}, Ptr{Void}, Clong, Clong,
							Cint, Clong, Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, d, t, p)
				 delete(unif01)
				 @printf "done "
				 return pvalue(sres)
		end
		function smarsa_GCD(gen::Generator, N::Integer, n::Integer, r::Integer, s::Integer)
				 unif01 = Unif01(gen, "")
				 sres = MarsaRes2(N)
				 ccall((:smarsa_GCD, libtestu01), Void,
						(Ptr{Void}, Ptr{Void}, Clong, Clong,
							Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, s)
				 delete(unif01)
				 return pvalue(sres)
		end
		function smarsa_CollisionOver(gen::Generator, N::Integer, n::Integer, r::Integer, d::Integer, t::Integer)
				 unif01 = Unif01(gen, "")
				 sres = MarsaRes()
				 ccall((:smarsa_CollisionOver, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
							Cint, Clong, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, d, t)
				 delete(unif01)
				 return pvalue(sres)
		end
		function smarsa_Savir2(gen::Generator, N::Integer, n::Integer, r::Integer, m::Integer, t::Integer)
				unif01 = Unif01(gen, "")
				sres = ResChi2(N)
				ccall((:smarsa_Savir2, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Clong, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, m, t)
				delete(unif01)
				return pvalue(sres)
		end
		function smarsa_SerialOver(gen::Generator, N::Integer, n::Integer, r::Integer, d::Integer, t::Integer)
				 unif01 = Unif01(gen, "")
				 sres = ResBasic()
				 ccall((:smarsa_SerialOver, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
							Cint, Clong, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, d, t)
				 delete(unif01)
				 return pvalue(sres)
		end
		function smarsa_MatrixRank(gen::Generator, N::Integer, n::Integer, r::Integer, s::Integer, L::Integer, k::Integer)
				@printf "\nTEST smarsa_MatrixRank \t"
				unif01 = Unif01(gen, "")
				sres = ResChi2(N)
				ccall((:smarsa_MatrixRank, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint, Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, s, L, k)
				delete(unif01)
				@printf "done "
				return pvalue(sres)
		end

		## sknuth
		function sknuth_Collision(gen::Generator, N::Integer, n::Integer, r::Integer, d::Integer, t::Integer)
				@printf "\nTEST sknuth_Collision \t"
				unif01 = Unif01(gen, "")
				sres = KnuthRes2()
				ccall((:sknuth_Collision, libtestu01), Void,
						(Ptr{Void}, Ptr{Void}, Clong, Clong,
								Cint, Clong, Cint),
						unif01.ptr, sres.ptr, N, n,
						r, d, t)
				delete(unif01)
				@printf "done "
				 return pvalue(sres)
		end
		function sknuth_CollisionPermut(gen::Generator, N::Integer, n::Integer, r::Integer, t::Integer)
				unif01 = Unif01(gen, "")
				sres = KnuthRes2()
				ccall((:sknuth_CollisionPermut, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, t)
				delete(unif01)
				return pvalue(sres)
		end
		function sknuth_CouponCollector(gen::Generator, N::Integer, n::Integer, r::Integer, d::Integer)
				@printf "\nTEST sknuth_CouponCollector \t"
				unif01 = Unif01(gen, "")
				sres = ResChi2(N)
				ccall((:sknuth_CouponCollector, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, d)
				delete(unif01)
				@printf "done "
				return pvalue(sres)
		end
		function sknuth_Gap(gen::Generator, N::Integer, n::Integer, r::Integer, Alpha::Real, Beta::Real)
				@printf "\nTEST sknuth_Gap \t"
				unif01 = Unif01(gen, "")
				sres = ResChi2(N)
				ccall((:sknuth_Gap, libtestu01), Void,
						(Ptr{Void}, Ptr{Void}, Clong, Clong,
								Cint, Float64, Float64),
						unif01.ptr, sres.ptr, N, n,
						r, Alpha, Beta)
				delete(unif01)
				@printf "done "
				return pvalue(sres)
		end
		function sknuth_MaxOft(gen::Generator, N::Integer, n::Integer, r::Integer, d::Integer, t::Integer)
				@printf "\nTEST sknuth_MaxOft \t"
				unif01 = Unif01(gen, "")
				sres = KnuthRes1(N)
				ccall((:sknuth_MaxOft, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, d, t)
				delete(unif01)
				@printf "done "
				return pvalue(sres)
		end
		function sknuth_Permutation(gen::Generator, N::Integer, n::Integer, r::Integer, t::Integer)
				unif01 = Unif01(gen, "")
				sres = ResChi2(N)
				ccall((:sknuth_Permutation, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, t)
				delete(unif01)
				return pvalue(sres)
		end
		function sknuth_Run(gen::Generator, N::Integer, n::Integer, r::Integer, up::Integer)
				unif01 = Unif01(gen, "")
				sres = ResChi2(N)
				ccall((:sknuth_Run, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, up)
				delete(unif01)
				return pvalue(sres)
		end
		function sknuth_SimpPoker(gen::Generator, N::Integer, n::Integer, r::Integer, d::Integer, k::Integer)
				@printf "\n sknuth_SimpPoker \t"
				unif01 = Unif01(gen, "")
				sres = ResChi2(N)
				ccall((:sknuth_SimpPoker, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, d, k)
				delete(unif01)
				@printf "done "
				return pvalue(sres)
		end

		## svaria
		function svaria_AppearanceSpacings(gen::Generator, N::Integer, Q::Integer, K::Integer, r::Integer, s::Integer, L::Integer)
				unif01 = Unif01(gen, "")
				sres = ResBasic()
				ccall((:svaria_AppearanceSpacings, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Clong, Cint, Cint, Cint),
						 unif01.ptr, sres.ptr, N, Q,
						 K, r, s, L)
				delete(unif01)
				return pvalue(sres)
		end
		function svaria_SampleProd(gen::Generator, N::Integer, n::Integer, r::Integer, t::Integer)
				unif01 = Unif01(gen, "")
				sres = ResBasic()
				ccall((:svaria_SampleProd, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, t)
				delete(unif01)
				return pvalue(sres)
		end
		function svaria_SampleMean(gen::Generator, N::Integer, n::Integer, r::Integer)
				unif01 = Unif01(gen, "")
				sres = ResBasic()
				ccall((:svaria_SampleMean, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r)
				delete(unif01)
				return pvalue(sres)
		end
		function svaria_SampleCorr(gen::Generator, N::Integer, n::Integer, r::Integer, k::Integer)
				unif01 = Unif01(gen, "")
				sres = ResBasic()
				ccall((:svaria_SampleCorr, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, k)
				delete(unif01)
				return pvalue(sres)
		end
		function svaria_SumCollector(gen::Generator, N::Integer, n::Integer, r::Integer, g::Float64)
				unif01 = Unif01(gen, "")
				sres = ResChi2(N)
				ccall((:svaria_SumCollector, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cdouble),
						 unif01.ptr, sres.ptr, N, n,
						 r, g)
				delete(unif01)
				return pvalue(sres)
		end
		function svaria_WeightDistrib(gen::Generator, N::Integer, n::Integer, r::Integer, k::Integer, alpha::Real, beta::Real)
				@printf "\nTEST svaria_WeightDistrib \t"	
				unif01 = Unif01(gen, "")
				sres = ResChi2(N)
				ccall((:svaria_WeightDistrib, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Clong, Float64, Float64),
						 unif01.ptr, sres.ptr, N, n,
						 r, k, alpha, beta)
				delete(unif01)
				@printf "done "
				return pvalue(sres)
		end

		## sstring
		function sstring_AutoCor(gen::Generator, N::Integer, n::Integer, r::Integer, s::Integer, d::Integer)
				unif01 = Unif01(gen, "")
				sres = ResBasic()
				ccall((:sstring_AutoCor, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, s, d)
				delete(unif01)
				return pvalue(sres)
		end
		function sstring_HammingCorr(gen::Generator, N::Integer, n::Integer, r::Integer, s::Integer, L::Integer)
				unif01 = Unif01(gen, "")
				sres = StringRes()
				ccall((:sstring_HammingCorr, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, s, L)
				delete(unif01)
				return pvalue(sres)
		end
		function sstring_HammingIndep(gen::Generator, N::Integer, n::Integer, r::Integer, s::Integer, L::Integer, d::Integer)
				@printf "\nTEST sstring_HammingIndep \t"
				unif01 = Unif01(gen, "")
				sres = StringRes()
				ccall((:sstring_HammingIndep, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint, Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, s, L, d)
				delete(unif01)
				@printf "done "
				return pvalue(sres)
		end
		function sstring_HammingWeight2(gen::Generator, N::Integer, n::Integer, r::Integer, s::Integer, L::Integer)
				unif01 = Unif01(gen, "")
				sres = ResBasic()
				ccall((:sstring_HammingWeight2, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
							Cint, Cint, Clong),
						 unif01.ptr, sres.ptr, N, n,
						 r, s, L)
				delete(unif01)
				return pvalue(sres)
		end
		function sstring_LongestHeadRun(gen::Generator, N::Integer, n::Integer, r::Integer, s::Integer, L::Integer)
				unif01 = Unif01(gen, "")
				sres = StringRes2(N)
				ccall((:sstring_LongestHeadRun, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
							Cint, Cint, Clong),
						 unif01.ptr, sres.ptr, N, n,
						 r, s, L)
				delete(unif01)
				return pvalue(sres)
		end
		function sstring_PeriodsInStrings(gen::Generator, N::Integer, n::Integer, r::Integer, s::Integer)
				unif01 = Unif01(gen, "")
				sres = ResChi2(N)
				ccall((:sstring_PeriodsInStrings, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, s)
				delete(unif01)
				return pvalue(sres)
		end
		function sstring_Run(gen::Generator, N::Integer, n::Integer, r::Integer, s::Integer)
				unif01 = Unif01(gen, "")
				sres = StringRes3(N)
				ccall((:sstring_Run, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
								 Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, s)
				delete(unif01)
				return pvalue(sres)
		end


		## swalk
		function swalk_RandomWalk1(gen::Generator, N::Integer, n::Integer, r::Integer, s::Integer, L0::Integer, L1::Integer)
				@printf "\nTEST swalk_RandomWalk1 \t"
				unif01 = Unif01(gen, "")
				sres = WalkRes(N)
				ccall((:swalk_RandomWalk1, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
							Cint, Cint, Clong, Clong),
						 unif01.ptr, sres.ptr, N, n,
						 r, s, L0, L1)
				delete(unif01)
				@printf	"done"
				return pvalue(sres)
		end

		## snpair
		function snpair_ClosePairs(gen::Generator, N::Integer, n::Integer, r::Integer, t::Integer, p::Integer, m::Integer)
				unif01 = Unif01(gen, "")
				sres = NpairRes(N)
				ccall((:snpair_ClosePairs, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
							Cint, Cint, Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, t, p, m)
				delete(unif01)
				return pvalue(sres)
		end

		## scomp
		function scomp_LempelZiv(gen::Generator, N::Integer, k::Integer, r::Integer, s::Integer)
				unif01 = Unif01(gen, "")
				sres = ResBasic()
				ccall((:scomp_LempelZiv, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Cint,
							Cint, Cint),
						 unif01.ptr, sres.ptr, N, k,
						 r, s)
				delete(unif01)
				return pvalue(sres)
		end
		function scomp_LinearComp(gen::Generator, N::Integer, n::Integer, r::Integer, s::Integer)
				unif01 = Unif01(gen, "")
				sres = CompRes(N)
				ccall((:scomp_LinearComp, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Clong,
							Cint, Cint),
						 unif01.ptr, sres.ptr, N, n,
						 r, s)
				delete(unif01)
				return pvalue(sres)
		end

		# sspectral
		function sspectral_Fourier3(gen::Generator, N::Integer, k::Integer, r::Integer, s::Integer)
				unif01 = Unif01(gen, "")
				sres = SpectralRes()
				ccall((:sspectral_Fourier3, libtestu01), Void,
						 (Ptr{Void}, Ptr{Void}, Clong, Cint,
							Cint, Cint),
						 unif01.ptr, sres.ptr, N, k,
						 r, s)
				delete(unif01)
				return pvalue(sres)
		end

		##################
		# Test Batteries #
		##################
		for (snm, fnm) in ((:SmallCrush, :smallcrushTestU01), (:Crush, :crushTestU01), (:BigCrush, :bigcrushTestU01), (:pseudoDIEHARD, :diehardTestU01), (:FIPS_140_2, :fips_140_2TestU01))
				@eval begin
						function $(fnm)(f::Generator, fname::ByteString)
								@printf "Unif01 call <---\n"
								unif01 = Unif01(f, fname)
								ccall(($(string("bbattery_", snm)), libtestu01), Void, (Ptr{Void},), unif01.ptr)
								delete(unif01)
						end
						$(fnm)(f::Generator) = $(fnm)(f::Generator, "")
				end
		end
		 function smallcrushJulia(f::Generator)
				 # @everywhere f = ()->g()   g->smarsa_BirthdaySpacings(g, 1, 5000000, 0, 1073741824, 2, 1),
				 testnames = [g->smarsa_BirthdaySpacings(g, 1, 5000000, 0, 1073741824, 2, 1),
											g->sknuth_Collision(g, 1, 5000000, 0, 65536, 2),
											g->sknuth_Gap(g, 1, 200000, 22, 0.0, .00390625),
											g->sknuth_SimpPoker(g, 1, 400000, 24, 64, 64),
											g->sknuth_CouponCollector(g, 1, 500000, 26, 16),
											g->sknuth_MaxOft(g, 1, 2000000, 0, 100000, 6),
											g->svaria_WeightDistrib(g, 1, 200000, 27, 256, 0.0, 0.125),
											g->smarsa_MatrixRank(g, 1, 20000, 20, 10, 60, 60),
											g->sstring_HammingIndep(g, 1, 500000, 20, 10, 300, 0)[:Mean],
											g->swalk_RandomWalk1(g, 1, 1000000, 0, 30, 150, 150)]
				 return pmap(t->t(f), testnames)
		 end
		 function bigcrushJulia(f::Generator)
				 testnames = [g->smarsa_SerialOver(g, 1, 10^9, 0, 2^8, 3)[:Mean],
											g->smarsa_SerialOver(g, 1, 10^9, 22, 2^8, 3)[:Mean],
											g->smarsa_CollisionOver(g, 30, 2*10^7, 0, 2^21, 2),
											g->smarsa_CollisionOver(g, 30, 2*10^7, 9, 2^21, 2),
											g->smarsa_CollisionOver(g, 30, 2*10^7, 0, 2^14, 3),
											g->smarsa_CollisionOver(g, 30, 2*10^7, 16, 2^14, 3),
											g->smarsa_CollisionOver(g, 30, 2*10^7, 0, 64, 7),
											g->smarsa_CollisionOver(g, 30, 2*10^7, 24, 64, 7),
											g->smarsa_CollisionOver(g, 30, 2*10^7, 0, 8, 14),
											g->smarsa_CollisionOver(g, 30, 2*10^7, 27, 8, 14),
											g->smarsa_CollisionOver(g, 30, 2*10^7, 0, 4, 21),
											g->smarsa_CollisionOver(g, 30, 2*10^7, 28, 4, 21),
											g->smarsa_BirthdaySpacings(g, 100, 10^7, 0, 2^31, 2, 1),
											g->smarsa_BirthdaySpacings(g, 20, 2*10^7, 0, 2^21, 3, 1),
											g->smarsa_BirthdaySpacings(g, 20, 3*10^7, 14, 2^16, 4, 1),
											g->smarsa_BirthdaySpacings(g, 20, 2*10^7, 0, 2^9, 7, 1),
											g->smarsa_BirthdaySpacings(g, 20, 2*10^7, 7, 2^9, 7, 1),
											g->smarsa_BirthdaySpacings(g, 20, 3*10^7, 14, 2^8, 8, 1),
											g->smarsa_BirthdaySpacings(g, 20, 3*10^7, 22, 2^8, 8, 1),
											g->smarsa_BirthdaySpacings(g, 20, 3*10^7, 0, 2^4, 16, 1),
											g->smarsa_BirthdaySpacings(g, 20, 3*10^7, 26, 2^4, 16, 1),
											g->snpair_ClosePairs(g, 30, 6*10^6, 0, 3, 0, 30),
											g->snpair_ClosePairs(g, 20, 4*10^6, 0, 5, 0, 30),
											g->snpair_ClosePairs(g, 10, 3*10^6, 0, 9, 0, 30),
											g->snpair_ClosePairs(g, 5, 2*10^6, 0, 16, 0, 30),
											g->sknuth_SimpPoker(g, 1, 4*10^8, 0, 8, 8),
											g->sknuth_SimpPoker(g, 1, 4*10^8, 27, 8, 8),
											g->sknuth_SimpPoker(g, 1, 10^8, 0, 32, 32),
											g->sknuth_SimpPoker(g, 1, 10^8, 25, 32, 32),
											g->sknuth_CouponCollector(g, 1, 2*10^8, 0, 8),
											g->sknuth_CouponCollector(g, 1, 2*10^8, 10, 8),
											g->sknuth_CouponCollector(g, 1, 2*10^8, 20, 8),
											g->sknuth_CouponCollector(g, 1, 2*10^8, 27, 8),
											g->sknuth_Gap(g, 1, 5*10^8, 0, 0.0, 1/16),
											g->sknuth_Gap(g, 1, 3*10^8, 25, 0.0, 1/32),
											g->sknuth_Gap(g, 1, 10^8, 0, 0.0, 1/128),
											g->sknuth_Gap(g, 1, 10^7, 20, 0.0, 1/1024),
											g->sknuth_Run(g, 5, 10^9, 0, 0),
											g->sknuth_Run(g, 5, 10^9, 15, 1),
											g->sknuth_Permutation(g, 1, 10^9, 0, 3),
											g->sknuth_Permutation(g, 1, 10^9, 0, 5),
											g->sknuth_Permutation(g, 1, 5*10^8, 0, 7),
											g->sknuth_Permutation(g, 1, 5*10^8, 10, 10),
											g->sknuth_CollisionPermut(g, 20, 2*10^7, 0, 14),
											g->sknuth_CollisionPermut(g, 20, 2*10^7, 10, 14),
											g->sknuth_MaxOft(g, 40, 10^7, 0, 10^5, 8),
											g->sknuth_MaxOft(g, 30, 10^7, 0, 10^5, 16),
											g->sknuth_MaxOft(g, 20, 10^7, 0, 10^5, 24),
											g->sknuth_MaxOft(g, 20, 10^7, 0, 10^5, 32),
											g->svaria_SampleProd(g, 40, 10^7, 0, 8)[:AD],
											g->svaria_SampleProd(g, 20, 10^7, 0, 16)[:AD],
											g->svaria_SampleProd(g, 20, 10^7, 0, 24)[:AD],
											g->svaria_SampleMean(g, 2*10^7, 30, 0)[:AD],
											g->svaria_SampleMean(g, 2*10^7, 30, 10)[:AD],
											g->svaria_SampleCorr(g, 1, 2*10^9, 0, 1)[:Mean],
											g->svaria_SampleCorr(g, 1, 2*10^9, 0, 2)[:Mean],
											g->svaria_AppearanceSpacings(g, 1, 10^7, 10^9, 0, 3, 15)[:Mean],
											g->svaria_AppearanceSpacings(g, 1, 10^7, 10^9, 27, 3, 15)[:Mean],
											g->svaria_WeightDistrib(g, 1, 2*10^7, 0, 256, 0.0, 1/4),
											g->svaria_WeightDistrib(g, 1, 2*10^7, 20, 256, 0.0, 1/4),
											g->svaria_WeightDistrib(g, 1, 2*10^7, 28, 256, 0.0, 1/4),
											g->svaria_WeightDistrib(g, 1, 2*10^7, 0, 256, 0.0, 1/16),
											g->svaria_WeightDistrib(g, 1, 2*10^7, 10, 256, 0.0, 1/16),
											g->svaria_WeightDistrib(g, 1, 2*10^7, 26, 256, 0.0, 1/16),
											g->svaria_SumCollector(g, 1, 5*10^8, 0, 10.0),
											g->smarsa_MatrixRank(g, 10, 10^6, 0, 5, 30, 30),
											g->smarsa_MatrixRank(g, 10, 10^6, 25, 5, 30, 30),
											g->smarsa_MatrixRank(g, 1, 5000, 0, 4, 1000, 1000),
											g->smarsa_MatrixRank(g, 1, 5000, 26, 4, 1000, 1000),
											g->smarsa_MatrixRank(g, 1, 80, 15, 15, 5000, 5000),
											g->smarsa_MatrixRank(g, 1, 80, 0, 30, 5000, 5000),
											g->smarsa_Savir2(g, 10, 10^7, 10, 2^30, 30),
											g->smarsa_GCD(g, 10, 5*10^7, 0, 30),
											g->swalk_RandomWalk1(g, 1, 10^8, 0, 5, 50, 50),
											g->swalk_RandomWalk1(g, 1, 10^8, 25, 5, 50, 50),
											g->swalk_RandomWalk1(g, 1, 10^7, 0, 10, 1000, 1000),
											g->swalk_RandomWalk1(g, 1, 10^7, 20, 10, 1000, 1000),
											g->swalk_RandomWalk1(g, 1, 10^6, 0, 15, 1000, 1000),
											g->swalk_RandomWalk1(g, 1, 10^6, 15, 15, 10000, 10000),
											g->scomp_LinearComp(g, 1, 400000, 0, 1),
											g->scomp_LinearComp(g, 1, 400000, 29, 1),
											g->scomp_LempelZiv(g, 10, 27, 0, 30)[:Sum],
											g->scomp_LempelZiv(g, 10, 27, 15, 15)[:Sum],
											g->sspectral_Fourier3(g, 100000, 14, 0, 3),
											g->sspectral_Fourier3(g, 100000, 14, 27, 3),
											g->sstring_LongestHeadRun(g, 1, 1000, 0, 3, 10^7),
											g->sstring_LongestHeadRun(g, 1, 1000, 27, 3, 10^7),
											g->sstring_PeriodsInStrings(g, 10, 5*10^8, 0, 10),
											g->sstring_PeriodsInStrings(g, 10, 5*10^8, 20, 10),
											g->sstring_HammingWeight2(g, 10, 10^9, 0, 3, 10^6)[:Sum],
											g->sstring_HammingWeight2(g, 10, 10^9, 27, 3, 10^6)[:Sum],
											g->sstring_HammingCorr(g, 1, 10^9, 10, 10, 30)[:Mean],
											g->sstring_HammingCorr(g, 1, 10^8, 10, 10, 300)[:Mean],
											g->sstring_HammingCorr(g, 1, 10^8, 10, 10, 1200)[:Mean],
											g->sstring_HammingIndep(g, 10, 3*10^7, 0, 3, 30, 0)[:Sum],
											g->sstring_HammingIndep(g, 10, 3*10^7, 27, 3, 30, 0)[:Sum],
											g->sstring_HammingIndep(g, 1, 3*10^7, 0, 4, 300, 0)[:Mean],
											g->sstring_HammingIndep(g, 1, 3*10^7, 26, 4, 300, 0)[:Mean],
											g->sstring_HammingIndep(g, 1, 10^7, 0, 5, 1200, 0)[:Mean],
											g->sstring_HammingIndep(g, 1, 10^7, 25, 5, 1200, 0)[:Mean],
											g->sstring_Run(g, 1, 2*10^9, 0, 3),
											g->sstring_Run(g, 1, 2*10^9, 27, 3),
											g->sstring_AutoCor(g, 10, 10^9, 0, 3, 1)[:Sum],
											g->sstring_AutoCor(g, 10, 10^9, 0, 3, 3)[:Sum],
											g->sstring_AutoCor(g, 10, 10^9, 27, 3, 1)[:Sum],
											g->sstring_AutoCor(g, 10, 10^9, 27, 3, 3)[:Sum]]
				 return pmap(t->t(f), testnames)
		 end
end

