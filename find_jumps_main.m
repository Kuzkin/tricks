% find_jumps_main

MAX_TIME_GAP = 10;  	  % separate data around such time gaps (sec) 
MIN_JUMP_SIZE = 100;    	% min jump sequence size points
jj1 = [];

Anrm = vnorm(ACC);
N = length(Ta);
dt = diff(Ta);
jj = find(dt>MAX_TIME_GAP);
jgaps = [1,jj',N];
M = length(jgaps);

zz = zeros(M-1,1);
AAav = zz;
AAstd = zz;
AAmax = zz;
Njmp = zz;
ii1 = zz; ii2=zz;
k1 = 0;

ALL_JUMPS=[];
for k=1:(M-1)
  i1 = jgaps(k)+2; i2 = jgaps(k+1)-2;
  if (i2-i1)>MIN_JUMP_SIZE
    ii = i1:i2;
    A = ACC(ii,:);
    T = Ta(ii);
	k1 = k1+1;
    figure(2);plot(Anrm(ii)); title(['k=',num2str(k1)])
	pause
    [all_jumps,Njmp(k1),AAav(k1),AAstd(k1),AAmax(k1)] = FindJumpsFuzzyFn4(A,T);
    ii1(k1) = i1; ii2(k1) = i2;
    [k1,Njmp(k1),AAav(k1),AAstd(k1),AAmax(k1)]
    if Njmp(k1)>0
        for i=1:Njmp(k1)
            j1 = all_jumps(i,1)+i1-1;
            j2 = all_jumps(i,2)+i1-1;
            ALL_JUMPS = [ALL_JUMPS;[j1,j2,Ta(j2)-Ta(j1)]];
        end
    end
    pause;
  end	
end
  
