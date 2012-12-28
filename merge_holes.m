% merge_holes.m

% mergin local holes			
% inside each potential jump merge segments with S1>0
% start merging from the right - from the Shock side
% the resulted merge must have at least one point with S1=1
% if there is a point with S1==1 then H(.)=1, otherwise H(.) = 0

Jump_start2 = zeros(size(Jump_start));
Jump_av = zeros(size(Jump_start));
for k=1:Nshk
	if Jump_start(k)>0
		i2 = Jump_end(k);
		i1 = Jump_start(k);
		i = i1+1;
		jm = 1;
		JH = [];
		JH(jm,1) = i1;
		H(jm) = 0;        % default
		while (i<i2)
			if S1(i)>0 & S1(i+1)==0		% end of local hole
				JH(jm,2) = i;
			end
			if S1(i)==0 & S1(i+1)>0		% start of local hole
				jm = jm+1;
				JH(jm,1) = i+1;
				H(jm) = 0;   		% so far no S1=1 value
			end
			i = i+1;
			if S1(i)==1 H(jm)=1; end	
		end
		JH(jm,2) = i2;
		Nhl = jm;
% now k-th potential Jump has Nhl small holes
% with indexes JH()
% if local hole m has S1==1 then H(m)==1
%		
% merge local holes from right to left based on the resulted average Asmz
% 
		Jstart_new = 0;
		Hflag = 0;				% no small points inside yet
		Jump_av(k) = 0;
		i2 = JH(Nhl,2);			% right edge of the hole
% the right most local hole always gets in, compute its Aav and Sav
		i1 = JH(Nhl,1);
		sav_old = mean(S1(i1:i2)); % max(mean(S1(i1:i2)), 1);
		aav_old = mean(Asmz(i1:i2));
		if H(Nhl)==1 			% if S==1 is in this hole
			Hflag = 1;
			Jstart_new = i1;
		end	
% now move to the left and merge if possible
% stop at the first non mergable local hole
		n = Nhl;
		merge = 1;
		while merge & n>1
			n = n-1;			% walk from right to left
			i1old = i1;
            i1 = JH(n,1);		% start of first local hole
           
			sav_new = mean(S1(i1:i2));
			aav_new = mean(Asmz(i1:i2));
            
            large_acc = max(Anrm(i1:i1old));
            Nlarge_acc = length(find(Anrm(i1:i1old) > 1.5*norm(Aav)));
            
           	ds = sav_new/sav_old; 		% should not go down too much
			da = aav_new/aav_old;		% should not go up too much
            
			if aav_new < Aav_max & large_acc < 1.5*norm(Aav) & Nlarge_acc < 5 &....
				ds > 0.5 & da < 2		% allow merge
				Jstart_new = i1;
				if H(n)==1 Hflag =1; end;
			else
				merge = 0;				% cannot merge, stop
			end
		end
		if Hflag
			Jump_start2(k) = Jstart_new;
			Jump_av(k) = aav_old;
		else
			Jump_start2(k) = 0;
		end
	end
end	
%
% Now potential jumps are defined by Jump_start2:Jump_end	