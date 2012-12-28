% find_holes.m
% find holes between istart:iend
% from left to right
% m-th hole coordinates are JH(m,1):JH(m,2)

	function [JH,H] = find_holes(S1,istart,iend)

		i2 = iend;
		i1 = istart;
		N = length(S1);
	% first make sure that you dont start and dont end in a middle of a hole
		s = S1(i1);
		while s>0 & i1>1
			i1 = i1 -1;
			s = S1(i1);
		end
		s = S1(i2);
		while s>0 & i2<N
			i2 = i2+1;
			s = S1(i2);
		end	
% now move from i1 to i2 to find all the holes		
		i = i1;
		jm = 0;
		JH = []; H =[];
		while (i<=i2)
			if S1(i)>0 & S1(i+1)==0		% end of local hole
				JH(jm,2) = i;
			end
			if S1(i)==0 & S1(i+1)>0		% start of local hole
				jm = jm+1;
				JH(jm,1) = i+1;
				H(jm) = 0;   			% so far no S1=1 value
			end
			i = i+1;
			if S1(i)==1 H(jm)=1; end	
		end

%% plot for debugging
%		[n,m] = size(JH);
%		ii = i1:i2;
%		jj = [];
%		for k=1:n
%			jj = [jj,JH(k,1):JH(k,2)];
%		end
%		plot(ii,S1(ii),ii,S1(ii),'.b',jj,S1(jj),'ro');	
			