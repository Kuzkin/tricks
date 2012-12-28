%   find_shocks.m

% find all suspicious shocks 
% (points where acceleration is between minshock and maxshock)

jm = 0;   % shock index
flag = 0; % flag in or out of shock
for k=1:N
   s = mf_0_to_1(Asmz(k),minshock,maxshock); % if s>0 then the point belongs to shock  
   Sshk(k) = s;
   if flag == 0
      if s>0    % start shock
         jm = jm+1;  % shock count
         Shk_start(jm) = k;
         flag = 1;	% inside shock
         shkw = 1;	% measure shock width
      end
   else    				% already inside shock  
      if s==0  % end of shock
         flag=0;
         Shk_w(jm) = shkw;
      else
         shkw = shkw + 1;
      end
   end
end
if s>0 Shk_w(jm) = shkw; end;    % if shock is not closed.
Shk_end = Shk_start+Shk_w-1; % array containing ends of shocks
Nshk = jm; % number of shocks

% find how sharp is each shock by expanding down left and right
% while drop larger than a threshold dAshk_min
 
for k=1:Nshk				% for all potential shocks
% first find peak and store it in Shk_top()
	i0 = Shk_start(k);
	amax = Asmz(i0);
	for i=(Shk_start(k)+1):Shk_end(k)
		if Asmz(i) > amax
			amax = Asmz(i);
			i0 = i;
		end
	end		
	Shk_top(k) = i0;
% now go left to find end of peak
    i = i0;
    a1 = Asmz(i);
	peak = 1;
   while (peak & i>1)
      i = i-1;
      a2 = Asmz(i);
	  if (a2>a1-dAshk_min) & a2<(Aav+2*Astd) peak=0; end
	  a1 = a2;
   end
   i1 = i+1;
% find right end of peak  
	i = i0;
    a1 = Asmz(i);
	peak = 1;
   while (peak & i<N)
      i = i+1;
      a2 = Asmz(i);
	  if (a2>a1-dAshk_min) & a2<(Aav+2*Astd) peak=0; end
	  a1 = a2;
   end
   i2 = i-1;
   Shk1(k) = i1;    	% Shock left edge
   Shk2(k) = i2;		% Shock right edge
   	shockhight = mf_0_to_1(Asmz(Shk_top(k)),minshock,maxshock); % how high is the shock
	w = Shk_top(k) - Shk1(k) +1;
	shocksharp_lft = mf_1_to_0(w/freq, minsharppeak_lft, maxsharppeak_lft); % is the shock sharp
	w = Shk2(k) - Shk_top(k) +1;
	shocksharp_rgh = mf_1_to_0(w/freq, minsharppeak_rgh, maxsharppeak_rgh);
	Shk_mf(k) = min([shockhight,shocksharp_lft,shocksharp_rgh]);
end
% now each shock is described by the following parameters
% Shk_top 		- index of the top of the shock
% Sshk	 		- shock initial mf - for all 1:N points
% Shk_h			- eqivalent hight in units of Astd
% Shk1, Shk2 	- points where shock starts and ends
% Shk_start, Shk_w, Shk_end - start, width, end of Sshk
% Shk_mf		- mf for the total shock - one number per shock
% 
