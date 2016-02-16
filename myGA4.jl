module myGA4
	using Base.AbstractRNG
	
	export SuperJuice, AbstractRNG
	
	abstract AbstractRNG
	
	type SuperJuice <:AbstractRNG
		function SuperJuice()
			n = nworkers()
			y = cell(n,2)
			y = pmap(long_computation3,32768,n,100,2,false,workers())
			ind = findmin(y[:,2])[2]
			return y[ind,1]
		end
		function SuperJuice(N,nc)
			n = nworkers()
			y = cell(n,2)
			y = pmap(long_computation2,1024,n,N,nc,false,workers())
			ind = findmin(y[:,2])[2]
                        return y[ind,1]
		end
			
	end

	function long_computation3(idx, N,kstop,nc,store2file)
		x = Array(Float64,N)
		filename = "long_run3_th"
		sigma0  = 0.485
		
		xc = Array(Float64,(N,nc))
		xc_score = ones(Float64,nc)
		ind_c = Array(Int64,(N,nc))
		
		xp =  Array(Float64,(N,nc)) #Array(Float64,N) 
		x_score::Float64 =0.0
		xp_score =ones(Float64,nc) #::Float64 =0.0
	
		x = EHEfast_0to1(rand(N),N)
		x_score = evaluate(x,N)
		
		xp_score = xp_score.*x_score;
		xc_score = xc_score.*x_score;
		
		if store2file
			csvfile = open(filename*string(myid())*".csv","w")
			csvfile2 = open(filename*string(myid())*"_extended"*".csv","w")
			score_tuple = tuple([x_score; xc_score; xp_score]...)
			write(csvfile,join( Any[0; score_tuple...], ","),"\n")
			write(csvfile2,join( x, ","),"\n")
		end

		#if (myid()==2)
		#	println("k=$(0) \t\t score = $(x_score)")
		#end

		k=1;
		while k<kstop
			score_tuple = tuple([x_score; xc_score; xp_score]...)
			while all(score_tuple[1].<= [score_tuple[2:end]...]) && k<kstop 
				for kk=1:nc
					ind_c[:,kk] = doOrder1X0_2(N)
					xc[:,kk] = copy(x[ind_c[:,kk]])
					xc[:,kk] = doShuffleMutation(xc[:,kk],N,0.1)
					xc_score[kk] = evaluate(xc[:,kk],N)
					
					xp[:,kk] = EHEfast_0to1(rand(N),N)
					xp_score[kk] = evaluate(xp[:,kk],N)
				end
				score_tuple = tuple([x_score; xc_score; xp_score]...)
				k +=1;
			end
			opt = findmin(score_tuple)[2]
			if opt>1 && opt<(nc+2)
				x = copy(xc[:,opt-1] )
				x_score = xc_score[opt-1]
			end
			if opt >= (nc+2)
				x = copy( xp[:,(opt - (nc+1)) ] )
				x_score = xp_score[(opt -(nc+1))]
			end
			
			
			if store2file
				write(csvfile,join( Any[k; score_tuple...], ","),"\n")
				write(csvfile2,join( x, ","),"\n")
			end
			
			#if (myid()==2)
			#	println("k = $(k) \t\t score = $(x_score)")
			#end

		end
		
		if store2file
			close(csvfile)
			close(csvfile2)
		end	
		
		return x, x_score
	end
			



	function long_computation1(idx, N,kstop,nc,store2file)
		x = Array(Float64,N)
		filename = "long_run3_th"
		
		xc = Array(Float64,(N,nc))
		xc_score = ones(Float64,nc)
		ind_c = Array(Int64,(N,nc))
		
		xp =  Array(Float64,(N,nc)) #Array(Float64,N) 
		x_score::Float64 =0.0
		xp_score =ones(Float64,nc) #::Float64 =0.0
	
		x = EHEfast_0to1(rand(N),N)
		x_score = evaluate(x,N)
		
		xp_score = xp_score.*x_score;
		xc_score = xc_score.*x_score;
		
		if store2file
			csvfile = open(filename*string(myid())*".csv","w")
			csvfile2 = open(filename*string(myid())*"_extended"*".csv","w")
			score_tuple = tuple([x_score; xc_score; xp_score]...)
			write(csvfile,join( score_tuple, ","),"\n")
			write(csvfile2,join( x, ","),"\n")
		end

		if (myid()==2)
			println("k=$(0) \t\t score = $(x_score)")
		end

		
		for k =1:kstop
			score_tuple = tuple([x_score; xc_score; xp_score]...)
			while all(score_tuple[1].<= [score_tuple[2:end]...]) 
				for kk=1:nc
					ind_c[:,kk] = doOrder1X0_2(N)
					xc[:,kk] = copy(x[ind_c[:,kk]])
					xc[:,kk] = doShuffleMutation(xc[:,kk],N,0.1)
					xc_score[kk] = evaluate(xc[:,kk],N)
					
					xp[:,kk] = EHEfast_0to1(rand(N),N)
					xp_score[kk] = evaluate(xp[:,kk],N)
				end
				score_tuple = tuple([x_score; xc_score; xp_score]...)
			end
			opt = findmin(score_tuple)[2]
			if store2file
				write(csvfile,join( score_tuple, ","),"\n")
			end	
			if opt>1 && opt<(nc+2)
				x = copy(xc[:,opt-1] )
				x_score = xc_score[opt-1]
			end
			if opt >= (nc+2)
				x = copy( xp[:,(opt - (nc+1)) ] )
				x_score = xp_score[(opt -(nc+1))]
			end
			if (myid()==2)
				println("k=$(k) \t\t score = $(x_score)")
			end
		end
		
		if store2file
			close(csvfile)
			
		end	
	
		return x, x_score
	end
			
	function long_computation0(idx, N,kstop,nc,store2file)
		filename = "long_run2_th"
		x = Array(Float64,N)
		
		xc = Array(Float64,(N,nc))
		xc_score = zeros(Float64,nc)
		ind_c = Array(Int64,(N,nc))
		
		xp = Array(Float64,N) 
		x_score::Float64 =0.0
		xp_score::Float64 =0.0
	
		x = EHEfast_0to1(rand(N),N)
		x_score = evaluate(x,N)

		
		if store2file
			csvfile = open(filename*string(myid())*".csv","w")
			csvfile2 = open(filename*string(myid())*"_extended"*".csv","w")
			score_tuple = tuple([x_score; xc_score; xp_score]...)
			write(csvfile,join( score_tuple, ","),"\n")
			write(csvfile2,join( x, ","),"\n")
		end
		if (myid()==2)
			println("k=$(0) \t\t score = $(x_score)")
		end

		for k =1:kstop
			for kk=1:nc
				ind_c[:,kk] = doOrder1X0_2(N)
				xc[:,kk] = copy(x[ind_c[:,kk]])
				xc[:,kk] = doShuffleMutation(xc[:,kk],N,0.1)
				xc_score[kk] = evaluate(xc[:,kk],N)
			end
			
			xp = EHEfast_0to1(rand(N),N)
			xp_score = evaluate(xp,N)
			
			score_tuple = tuple([x_score; xc_score; xp_score]...)
			opt = findmin(score_tuple)[2]
				
			if opt>1 && opt<(nc+2)
				x = copy(xc[:,opt-1] )
				x_score = xc_score[opt-1]
			end
			
			if opt == (nc+2)
				x = copy(xp)
				x_score = xp_score
			end
			
			if store2file
				write(csvfile,join( score_tuple, ","),"\n")
				if (k==10) || (k==100) || (k == 500) || (k == 1000) || (k== 5000) || (mod(k,10000)==0) 
					write(csvfile2,join( x, ","),"\n")
				end
			end
			if (myid()==2)
				println("k = $(k) \t\t score = $(x_score)")
			end
		end

		if store2file
			close(csvfile)
			close(csvfile2)
		end	
	
		return x, x_score
	end

	function EHEfast_0to1(xin,N)
		ind = sortperm(xin)
		ind2 = zeros(N)
		[ind2[ind[k]] = k for k=1:N]
		binsize = 2
		ratio = convert(Float64, binsize)/convert(Float64, N)
		offset = floor((ind2-1)/binsize)*ratio
		xin = mod(xin,ratio) + offset
		return xin
	end

	function evaluate(S,N)
		
		s = [S-mean(S); zeros(Float64,N)]
		temp = fft(s);		
                psd = (1.0/(Float64(2*N-1)*var(S) ) ).* (  temp.*conj( temp ) )
		return std(psd)
	end

	function doOrder1X0_2(N)
		p1_ind = shuffle([1:N;])
		p2_ind = shuffle([1:N;])
		childVec_ind = zeros(Int64,N)
		ini_Pt = rand( [1:(N - 1 );]);
		end_Pt = rand( [(ini_Pt + 1):N;]);

		A1 = zeros(Int64, end_Pt-ini_Pt + 1 )
		A1 = p1_ind[ini_Pt:end_Pt]
		A0 = setdiff(p2_ind, A1 )[1:(ini_Pt-1)]
		A2 = setdiff(p2_ind, [A1;A0] )
		childVec_ind = [A0; A1; A2]
		return childVec_ind
	end
	
	function doShuffleMutation(child1, N, prob_mutation)
		if rand() <= prob_mutation
			child1 = shuffle(child1)
		end
		return child1
	end
	#(long_computation0,1024,n,N,nc,true,workers())
	#(p, f, idx, N, kstop, nc, save2file)
	function pmap(f, N, M, kstop, nc, save2file, myprocs=workers())
		np = nprocs()  # determine the number of processes available
		results = cell(M,2)
		i = 1
		# function to produce the next work item from the queue.
		# in this case it's just an index.
		nextidx() = (idx=i; i+=1; idx)
		@sync begin
		   for p in myprocs
		       if p != myid() || np == 1
		           @async begin
		               while true
		                   idx = nextidx()
		                   if idx > M
		                       break
		                   end
		                   results[idx,1],results[idx,2]  = remotecall_fetch(p, f, idx, N, kstop, nc, save2file)
		               end
		           end
		       end
		   end
		end
		return results
	end

end
	

