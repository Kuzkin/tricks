 % ------ Main variables ------
% N - number of points in the whole run
% ACC - raw accelerometer readings
% A - cutted ACC
% Anrm - acceleration norm
% Asmz - smoothed acceleration norm
% freq - frequency of data collection
% S - array containing current values of jump membership function
% Sshk - array containing values of shock membership function
% T - GPS time
% Ta - sensor time
% V - GPS velocity (inplane components)

function [ALL_JUMPS,Njumps,Aav,Astd,Amax] = FindJumpFuzzyFn4(ACC,Ta)

N = length(ACC);

% historically in some mat files time is Tacc and in some Ta
% use Ta in the code
if exist('Tacc')
   if length(Tacc)==length(ACC)
      Ta = Tacc;
   end
end

% determine frequency
if (exist('Ta'))
   %freq =  1/median(diff(Ta));    % median is too expensive
   dt=diff(Ta(1:min(500,length(Ta))));
   jj = find(dt<0.5);
   freq = 1/mean(dt(jj));
else
    freq = 66;
end;

% remove deltafunctions, smooth data
h = 1/freq;
A = ACC(1:N,:);
N = length(A);
ACC0=A;
% if freq is low - then no smoothing -> asmz=1
% max smoothing is 0.5
if freq<=25  asmz = 1; end
if freq>=50	 asmz = 0.5; end
if freq<50 & freq>25 asmz = 1 -0.02*(freq-25); end
%asmz = 0.5+ 0.5 * mf_0_to_1(freq,25,50);
Anrm = vnorm(A);

%fnd_delta; 
rmv_delta; 						% removing delta-like peaks
Asmz = smz_exp(Anrm,asmz); % smooth data using exponential smoothing

% compute average and std of active points for scaling
if (~exist('ii_active')) 
    ii_active=1:N; 
end;
Aav = mean(Asmz(ii_active));
Astd = std(Asmz(ii_active));
Amax = max(Asmz(ii_active));

% ------------------------ MAGIC CONSTANTS for Fuzzy algorithm -----------------------
% maxsusp = 7;
maxshock = Aav + 6*Astd; % if ACC>maxshock then it is a shock for sure.
if (maxshock < 40) 
    maxshock = 40; 
end;

if (Aav + 4*Astd > Amax - 2*Astd)
    minshock = Aav + 2*Astd;
else
    minshock = Aav + 4*Astd;
end;

minsharppeak_lft = 3*0.020; 		% if <time in ms then mf=1
maxsharppeak_lft = 15*0.020;		% if > time in ms then mf=0
minsharppeak_rgh = 5*0.020; 		% if <time in ms then mf=1
maxsharppeak_rgh = 30*0.020;		% if > time in ms then mf=0
maxshockdist 	 = 15*0.020; 		% if time in ms between hole and jump > X then mf=0
maxjump2dist	 = 0.3;				% dist to jump must be < jump width
maxbridge2shock  = 0.5;				% total ACC in shock must be > than in bridge
minjumpsize = 0.15;   	% min duration of jump in sec
maxjumpsize = 2.00;  	% if duration > maxjumpsize (sec) then it is a jump for sure.
maxsitdowns = 5; % max sitdowns in the vicinity of the set. If number of sets is > maxsitdowns then it is not a jump
maxdistsitdowns = 500; % number of points to the right and to the left from the SET that may contain sitdowns
minjumptositratio = 2; % real jump should be at least minjumptositratio times longer than nearby sitdowns
maxjumptositratio = 3; % if jump is maxjumptositratio times longer than nearby sitdowns then it is a jump for sure
smallacc_smz = 4;      % the value of smoothed acceleration that is consedered as small
maxsusp_smz = Aav - 1.0*Astd;
Tjump_max = 2.0;  % max jump 2 sec
dAshk_min = 5;		% min change in A to be shock
Aav_max = Aav - 0.5*Astd;	% max allowed mean(Asmz) over jump; used for hole merging
% ------------------------ eo MAGIÑ CONSTANTS for Fuzzy algorithm -----------------------


% compute mf (S) for all points ; mf = 1 for Asmz<small_acc_smz
S = zeros(N,1);   % matlab is more efficient if array is created
for k=1:N
   S(k) = mf_1_to_0(Asmz(k), smallacc_smz, maxsusp_smz);
end

S0=S;     % initial value of mf

% find all suspicious shocks
Sshk = zeros(N,1);    		% shock mf 
find_shocks;

% now each shock is described by the following parameters
% Shk_top 		- index of the top of the shock
% Sshk	 		- shock initial mf - for all 1:N points
% Shk_h			- eqivalent hight in units of Astd
% Shk1, Shk2 	- points where shock starts and ends
% Shk_start, Shk_w, Shk_end - start, width, end of Sshk
% Shk_mf		- mf for the total shock - one number per shock
% 

% for each shock find holes left of it
% search only Tjump_max sec before shock
% however, stop before another shock

ndist = fix(Tjump_max*freq) + 1;		% num of points for Tjump_max
S1 = zeros(N,1);

find_possible_jumps;

% now suspected jumps are at S1>0 
% they start at Jump_start(k) and end at Jump_end(k)
% with shock at Shk_start(k)

% % show all the holes and shocks
% ii=ii_active;
% jj=find(Jump_start>0);  kk=Jump_start(jj);  kk2=Jump_end(jj);
% mm = Shk_start;
% figure(1)
% plot(ii,Asmz(ii),ii,5*S1(ii),'r',kk,Asmz(kk),'.g',...
%                 kk2,Asmz(kk2),'.m',mm,Asmz(mm),'r.')

% find sitdowns using S1 function
find_sitdowns;
           
% mergin local holes			
% inside each potential jump merge segments with S1>0
% start merging from the right - from the Shock side
% the resulted merge must have at least one point with S1=1

merge_holes;

% now potential jumps are defined by Jump_start2:Jump_end
% average value of Asmz over jump is Jump_av()
% mf function is still in S1
	
% plot each shock - hole pair to evaluate merging results
% green 		- initial potential Jump points
% green dots	- potential Jump point after merging
% red circle 	- shock_top
% magenta dots 	- shock peak points
% also jump mf S1

% figure(2);
%plt_shocks;

% up to this point we worked with Shocks and Jumps separately.
% now we need to decide which pair of jump-shock to keep
% parameter to consider:
% shock mf, max, sharpness, and equiv height
% hole mf, width, Aav, and how far it is from shock

Smerged = S1;

% The next step is to check all SETs and decide whether the SET can be
% considered as a jump or not. 
W = [];
sitdown_mf = [];
for k=1:Nshk
        sit_difference_mf(k) = 0;
        sit_width_mf(k) = 0;
		if Jump_start2(k) == 0
			Hole_mf(k) = 0;
            Nsit_mf(k) = 0;
            Nsits(k) = 0;
            angleLR(k) = 0;
			Hole_all(k,:) = [0,0,0];
            
            average_shock_mf(k)=  0;
            good_shock_mf(k)   =  0;
            average_hole_mf(k) =  0;
            good_hole_mf(k)    =  0;  
            Hole_and_Shock_mf(k) = 0;
        else	
			i1 = Jump_start2(k);
			i2 = Jump_end(k);
			d = Shk1(k) - Jump_end(k) - 1;
			shockdist = mf_1_to_0(d/freq,0,maxshockdist); %  the shock is far from the SET           
            w = i2 -i1 + 1;
            W(k) = Ta(i2) - Ta(i1);
	% compare jump width and dist from jump to shock - must be dist<<width		
			jump2distw = mf_1_to_0(d/w,0,maxjump2dist);
	% accel between jump and shock must be smaller than shock
			ashock = sum(Asmz(Shk1(k):Shk2(k)));
			j1 = Jump_end(k)+1;
			j2 = Shk1(k)-1;
			if j2<j1
				abridge = 0;
				bridge2shock = 1;
			else
				abridge = sum(Asmz(j1:j2));
				bridge2shock = mf_1_to_0(abridge/ashock,0,maxbridge2shock);
            end;
	% shock bridge must be small OR narrow
			shockdist_mf = min([shockdist,max([jump2distw,bridge2shock])]);
			w = mean(Asmz(i1:i2));
			
            jumpaccav = 1; % mf_accav(w,smallacc_smz,maxacc_smz);
			jumpsav = mean(S1(i1:i2));
			Hole_all(k,:) = [shockdist_mf,jumpaccav,jumpsav];
			Hole_mf(k) = min(Hole_all(k,:));
            
% --------------------------------------------------------------------------
        
            average_shock_mf(k)=  mf_1_to_0(Shk_mf(k), 0, 0.5);
            good_shock_mf(k)   =  mf_0_to_1(Shk_mf(k), 0.25, 0.5);
            average_hole_mf(k) =  mf_1_to_0(W(k), 0.25, 0.75);
            good_hole_mf(k)    =  mf_0_to_1(W(k), 0.5, 1);  
            Hole_and_Shock_mf(k) = max([min([Shk_mf(k), Hole_mf(k)]),  min([average_shock_mf(k), good_hole_mf(k)]),... 
                                        min([good_shock_mf(k), average_hole_mf(k)]), min([good_shock_mf(k), good_hole_mf(k)])]);
                                    
            if (W(k)/(Jump_end(k) - Jump_start2(k) + 1) > 0.045)
                Hole_mf(k) = 0; Shk_mf(k) = 0;
                Hole_and_Shock_mf(k) = 0;
            end; 
        end	
        
 % Check if there are a lot of sitdowns around the set.
 % If yes then it is not a jump
 % Calculate sitdowns around the hole.
 % Three sets are checked for sitdowns: maxdistsitdowns points to 
 % the left from the hole, maxdistsitdowns points to the right 
 % from the hole and 0.5*maxdistsitdowns points in both directions 
 
 calculate_sitdowns; 

% CURRENTLY NOT USED  
% %   determine relative depth of the hole
% %   if the hole is not deep in comparison with the surroundings 
% %   then it is not a jump 
% 
%  depth_mf(k) = 1;
%  if (istart > 200 & N - iend > 200) 
%        i01 = istart;
%        i02 = iend;
%        ii = i01:i02;
%        av = mean(Anrm(ii));
%        i1 = istart-100;
%        i2 = iend+100;
%        ss = sort(Anrm(ii));
%        i = fix(0.8*length(ii));
%        Ain(k) = ss(i);
%        ii_out = [(istart-200):(istart-1),(iend+1):(iend+200)];
%        Aout(k) = mean(Anrm(ii_out));
%        aa = sort(Anrm(ii_out));
%        i = fix(0.8*length(ii_out));
%        Aout2(k) = aa(i);
% 
%        if (10*Ain(k)/Aout(k) > 7) 
%            depth_mf(k) = 0;
%        end;
% end;
% end of CURRENTLY NOT USED  

 
% Calculate jump membership function
% Jump is Shock AND Hole AND sitdowns_mf
% ------------------------------------------------------------------------
 
sitdown_mf(k) = finalsitdowns_mf;
Jump_mf(k) = min([Hole_and_Shock_mf(k), sitdown_mf(k)]); % 

end;

show_jumps;


