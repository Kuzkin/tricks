% Calculate sitdowns around the hole.
% Array sitmiddles defined in function find_sitdowns is used.
% Three sets are checked for sitdowns: maxdistsitdowns points to 
% the left from the hole, maxdistsitdowns points to the right 
% from the hole and 0.5*maxdistsitdowns points in both directions. 
        
istart = Jump_start2(k); 
iend   = Jump_end(k); 

leftbound  = maxdistsitdowns;
rightbound = maxdistsitdowns;
if (iend ~=0 & istart ~= 0 & istart < N)
        ileft = istart;
        flag = 0;
        while (ileft > 1 &  flag == 0)
           if (abs(Ta(ileft) - Ta(ileft-1)) < 2 * h) 
               ileft = ileft-1;
           else flag = 1;
           end;
        end;
        leftbound = istart - ileft;

        iright = iend;
        flag = 0;
        while (iright <  N &  flag == 0) 
            if(abs(Ta(iright+1) - Ta(iright)) < 2 * h)
                iright = iright + 1;
            else flag = 1;
            end;
        end;
        rightbound = iright - iend;
end;
%-------------------------------------------------------------------   
% calculate number of sit downs N_sit around the jump

b1 = istart - min(maxdistsitdowns, leftbound);
b2 = iend + min(maxdistsitdowns, rightbound);
b3 = istart - 0.5*min(maxdistsitdowns, leftbound);
b4 = iend + 0.5*min(maxdistsitdowns, rightbound);
  
sitleft   = find(sitmidle > b1 & sitmidle < istart);
sitmid    = find(sitmidle > iend  & sitmidle < b2);
sitright1 = find(sitmidle > b3 & sitmidle < istart);
sitright2 = find(sitmidle > iend & sitmidle < b4);  

Nsits(k) = max([length(sitleft), length(sitright1)+length(sitright2), length(sitmid)]); 


mindist = 0;
centralsit = find(sitmidle > istart & sitmidle < iend);
    
w_jump = 0;
if (length(centralsit) > 0)
      w_jump = sitlength(centralsit(1));
end;

% calculate "number" of sit downs N_sit more accurately. The contribution
% of sit down to N_sit depends on its length. Sit downs are usually not
% longer than 0.2 sec.

sitl = [];
if (Nsits(k) >= 3)
        avdistbtjumps =  maxdistsitdowns/Nsits(k); 
        maxdistsitdowns1 = maxdistsitdowns + avdistbtjumps;
        b1 = max(1, istart - maxdistsitdowns1);
        b2 = min(N, iend + maxdistsitdowns1);
        b3 = max(1, istart - 0.5*maxdistsitdowns1);
        b4 = min(N, iend + 0.5*maxdistsitdowns1);

        sitleft2  = find(sitmidle > b1 & sitmidle < istart);  
        sitmid2   = find(sitmidle > iend  & sitmidle < b2);
        sitright21 = find((sitmidle > b3 & sitmidle < istart));
        sitright22 = find(sitmidle>iend & sitmidle < b4);
                   
        N_sitleft = 0;        N_sitmid = 0;        N_sitright = 0;

        kk = 1;
        for (ii=1:length(sitleft2))
              sitl(kk) = sitlength(sitleft2(ii));
              N_sitleft = N_sitleft +  mf_1_to_0(sitl(kk), 0.2, 0.5);
              kk = kk+1;
        end;
        for (ii=1:length(sitright21))
              sitl(kk) = sitlength(sitright21(ii));
              N_sitright = N_sitright + mf_1_to_0(sitl(kk), 0.2, 0.5);
              kk = kk+1;
        end;
        for (ii=1:length(sitright22))
              sitl(kk) = sitlength(sitright22(ii));
              N_sitright = N_sitright + mf_1_to_0(sitl(kk), 0.2, 0.5);
              kk = kk+1;
        end;
        for (ii=1:length(sitmid2))
              sitl(kk) = sitlength(sitmid2(ii));
              N_sitmid = N_sitmid + mf_1_to_0(sitl(kk), 0.2, 0.5);
              kk = kk+1;
        end;
        Nsits(k) = max([N_sitleft, N_sitright, N_sitmid]);
end;

Nsit_mf(k) =  mf_1_to_0(Nsits(k), -3.5, maxsitdowns);     

% check if the hole is the same as neighboring  sitdowns
% it is not the same if mean length of the sitdown is much smaller 
    if (sitl > 0)
        if (w_jump>0)
            sit_difference_mf(k) = mf_0_to_1(w_jump/mean(sitl), minjumptositratio, maxjumptositratio);
        else sit_difference_mf(k) = 1;
        end;
    end;
  
 finalsitdowns_mf = max(Nsit_mf(k), sit_difference_mf(k));    
 