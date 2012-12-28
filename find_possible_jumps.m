% find_possible_jumps.m

% for each shock find potential jump that is left from shock
% compute mf_suspicious for all points 2 sec left of shock
% if there are low enough measurements - then it is a potential jump
% then remove mf=0 at the ends of the interval and the rest is a potential jump for this shock

for k=1:Nshk				% for each shock
   i2 = Shk_start(k);		% end point of seach for holes - left edge of shock
   if k==1 
      ishk = 1;
   else
      ishk = Shk_end(k-1); 	% dont go farther than nearest left shock
   end
   i1 = max([1,i2-ndist,ishk]);  % i1:i2 potential jump
% make sure that jump has mf==1 (and that there is no shock inside)
	jmp = 0; 
   for i=i1:i2
	S1(i) = mf_1_to_0(Asmz(i), smallacc_smz, maxsusp_smz);
    if S1(i)==1 jmp = 1; end	% found low point
   end
   if jmp==0
      Jump_start(k) = 0;    % there is no jump here
      Jump_end(k) = 0;
   else    
% remove all points with mf=0 on left side
		i = i1;
		s = S1(i);
		while ( s==0 & i<i2)
			i = i+1;
			s = S1(i);
		end 
		Jump_start(k) = i;
% remove all points with mf=0 on right side
		i = i2;
		s = S1(i);
		while ( s==0 & i>i1)
			i = i-1;
			s = S1(i);
		end 
		Jump_end(k) = i;
   end    
end

% now suspected jumps are at S1>0 only
% they start at Jump_start(k) and end at Jump_end(k)
% with shock at Shk_start(k)
