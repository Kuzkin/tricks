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
   jj = find(dt<0.1);
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
if (maxshock<40) 
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
