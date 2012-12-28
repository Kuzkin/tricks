% show_jumps.m

TJMP_MIN = 300;   	% min jump size (ms) at mf=0
a0 = 1.5/TJMP_MIN;  	% 1.5 is a magic constant
Ain2Aout_MAX = 0.6;	% max allowed ratio of Anrm inside and outside of jump
Ain2Aout_ABS = 0.8; % ultimate max allowed ratio of Anrm inside and outside of jump
MF_MIN		= 0.1;	% min mf allowed for a jump
MIN_HOLE_VOL = 2000;
PASS_HOLE_VOL = 5500;
HOLE_M = 4;          % search hole border 4 points left and right
TRUE = 1;
FALSE = 0;

Njumps = 0;				% default init
ALL_JUMPS = [];

jj = find(Jump_mf>0 & Jump_start2>0 & Jump_end>0);
W = (Ta(Jump_end(jj)) - Ta(Jump_start2(jj)))*1000;
mm = find(W>0 & W<5000);  % avoid jumps between runs
W = W(mm);
jj = jj(mm);

% create cutout region boundaries in (W,mf) plane
y = -a0*(W'-TJMP_MIN);
kk = find(y<MF_MIN); y(kk) = MF_MIN;

kk1 = find(Jump_mf(jj)>y);    % jumps that passed
kk2 = find(Jump_mf(jj)<=y);	  % jumps that rejected

i_jump = jj(kk1);			  % jump index that passed (W,mf) region cutoff

n = length(i_jump);
if n==0 sprintf('No Jumps Found'), return; end

%% now all the jump candidates that pass (w,mf) region test are in i_jump() array
%% the last test on how "cup like" whole looks

% prepare (x,y) to plot region boarders on (w,mf) plane

x = 0:TJMP_MIN; 
y = -a0*(x-TJMP_MIN);
%w = 0:(fix(max(W))+1);
w = [0;W];
z = MF_MIN*ones(size(w));
zero = zeros(size(Anrm));

Ajumps = [];
Anew = []; Ain = [];
Av_out = [];
r_in_out = []; r_in_out2 = [];

A2out = [];

% for each jump that passed plot Anrm() around jump
% compute Ain/Aout and "volume"
% and show W vs mf

% put results into array ALL_JUMPS
% [start_indx, end_indx, width_ms]
ALL_JUMPS=[];

pass_flag = [];                 % last test result
n_jumps_good = 0;				% number of jumps that passed all tests
n_jumps_bad = 0;
k_jumps_good = [];
k_jumps_bad = [];



for k=1:n
   j = i_jump(k);				% index in jump arrays
   i01 = Jump_start2(j);		% jump start index
   i02 = Jump_end(j);			% jump end index
   mf = Jump_mf(j);				% jump mf
   m = kk1(k);					% index in jumps that passed (w,mf) region inclusion test
   ww = W(m);					% jump width in ms
   ii = i01:i02; 				% points inside jump
 % compare 80% value inside the jump with aver outside near jump  
 % we should check that the range is inside one active segment
 % but this is done in C++ so no need to complicate code here
   Ain(k) = mean(Anrm(ii));   % for debug only
   ss = sort(Anrm(ii));
   i = max([1,fix(0.8*length(ii))]);
   Anew(k) = ss(i);	   % 80% value
   i1_out = max([1,i01-200]);
   i2_out = min([N,i02+200]);
   ii_out = [i1_out:(i01-1),(i02+1):i2_out];   % points left and right from jump
   Av_out(k) = mean(Anrm(ii_out));	% aver outside
   r_in_out(k) = Anew(k)/Av_out(k);
   r_in_out2(k) = Ain(k)/Av_out(k);  % for debug only

   % now compute the "volume of water" that jump hole can hold
   % and do the last test
   hShk = Anrm(Shk_top(j));
   [vol(k),jh] = cmp_hole_volume2(Anrm,ww,hShk,i01,i02,HOLE_M);
  
   % -------------------   compute final Membership Function   ---------
   Hole_and_Shock_mf(j)
   shock = Shk_mf(j);
   hole = Hole_mf(j);
   sitdown = sitdown_mf(j);
   volume_mf = mf_0_to_1(vol(k), MIN_HOLE_VOL, PASS_HOLE_VOL);

   in_out = mf_1_to_0(r_in_out(k), Ain2Aout_MAX - 0.1, Ain2Aout_MAX + 0.1);
   
   if (vol(k) == 0) in_out = 0;  end;
   
   average_shock =  mf_1_to_0(shock, 0, 0.5);
   good_shock    =  mf_0_to_1(shock, 0.25, 0.5);
   average_width =  mf_1_to_0(0.001*ww, 0.25, 0.75);
   good_width    =  mf_0_to_1(0.001*ww, 0.5, 1);  
   
   % last test by hole shape
   if ((r_in_out(k) < Ain2Aout_MAX & vol(k) > MIN_HOLE_VOL )|... 
       (r_in_out(k) < Ain2Aout_ABS & vol(k) > PASS_HOLE_VOL ))
      shape = 1;
   else
	  shape = 0;
   end
   
   widht(k)    = ww;
   Final_mf(k) = max...
   ([...
       min([mf, shape]),...  % average width, average shock, but good shape
       min([good_shock, average_width, sitdown, in_out]),... % average jump with good shock, good hole, but not a sitdown
       min([average_shock, good_width, volume_mf]),... % long jump with bad shock, but good hole
       min([good_shock, good_width, volume_mf])... % long jump with good shock and hole  
   ]);
 
% check if the final mf is inside available region on (w,mf) plot
   ih = i01-HOLE_M+jh-1;
   if (Final_mf(k) > max(-a0*(ww - TJMP_MIN), MF_MIN))
      pass_flag(k) = TRUE;
	  n_jumps_good = n_jumps_good + 1;
	  k_jumps_good(n_jumps_good) = i_jump(k);
	  ALL_JUMPS(n_jumps_good,:) = [i01,i02,ww];
   else
	  pass_flag(k) = FALSE;
	  n_jumps_bad = n_jumps_bad + 1;
	  k_jumps_bad(n_jumps_bad) = i_jump(k);
   end

%-------------------------------------------------------------------

   % here compute all kind of stats for investigation
   % these stats are not used in jump filterming
   av = mean(Anrm(ii));					% aver value of Anrm inside jump
   avr = zero; avr(ii)=av;
   amx = max(Anrm(ii));
   aa = sort(Anrm(ii_out));
   i = max([1,fix(0.8*length(ii_out))]);
   A2out(k) = aa(i);
   
   % all computations are done, plot this jump
   
   % --------------------   Vertical and Magnetic  ------------------------------   
   
   if(exist('Mnormalized'))
     rotation_magnetic;
     magangle; 
	 MAG_DATA = 1;
   else
	 MAG_DATA=0;
   end
 % -------------------- eo Vertical and Magnetic  ------------------------------
 
% plot data around jump
 
  plt_jump_mag;
 
 % plot old and final mf vs. width for all jumps
  
  plt_jump_mf_w;
  
  pause
%%   Ajumps = [Ajumps;0;0;0;Anrm(ii)];		% array of all jumps put together
end

% end of computation
% there are n_jumps_good jumps that passed
% their stats are in ALL_JUMPS
   
Njumps = n_jumps_good;
sprintf('Total Found %d Jumps',Njumps)

% now plot all the jump on one (w,mf) plot

plt_all_jump_mf_w;

 

   
