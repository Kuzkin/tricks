% % filename   =  ('data\only_jumps_with_mag');
% % JUMP 1
% % Answer - 164 deg

% JUMP 2
% pstart = 800;
% pend   = 900;
% shock1 = 1250;
% shock2 = 1450;
% % Answer - 162

% % JUMP 3
% pstart = 1250;
% pend   = 1400;
% shock1 = 1686;
% shock2 = 1696;
% % Answer - 30

% % % JUMP 4
% pstart = 2100;
% pend   = 2200;
% shock1 = 2440;
% shock2 = 2470;
% % Answer - 35

% % % % JUMP 5
% pstart = 2600;
% pend   = 2700;
% shock1 = 2950;
% shock2 = 2970;
% % % % Answer - 30

% % % % % JUMP 6
% pstart = 3800;
% pend   = 3900;
% shock1 = 4190;
% shock2 = 4210;
% % % % % Answer - 14  

% pstart  - start of plateau
% pend    - end of plateau
% shock1  - start of shock or plateau after jump
% shock2  - end of shock or plateau after jump


% --------------------   plateau before jump and shock -----------------
pstart = i01_prev;
pend   = i01;

shock1 = Shk_start(jshk);
shock2 = Shk_end(jshk);

% ---------------  BEFORE JUMP ----------------------
% vertical before jump
Vertical1 = mean(A(pstart:pend,:))/norm(mean(A(pstart:pend,:)));  % vertical before jump
eM_before = [];
eA_before = [];
for (ii=pstart:pend)
    if (Anrm(ii) < 12.5 && Anrm(ii) > 7.5)
        eA_before = [eA_before;  A(ii,:)/norm(A(ii,:))];
        eM_before = [eM_before;  Mnormalized(ii,:)/norm(Mnormalized(ii,:))];
    end;
end;
Vertical1 = mean(eA_before); Vertical1 = Vertical1/norm(Vertical1);
eM_mean = mean(eM_before); eM_mean = eM_mean/norm(eM_mean);

anglebefore_ = acos(eM_mean*Vertical1')*180/3.14;

% plane projection of magnetic before jump
magplane = [];
for (ii=pstart:pend)
    if (Anrm(ii) < 12.5 && Anrm(ii) > 7.5)
        magplane = [magplane; Mnormalized(ii,:) -  (Mnormalized(ii,:)*Vertical1')*Vertical1];
    end;
end;
meanmagbefore =  mean(magplane)/norm(mean(magplane));

% ---------------  AFTER JUMP ----------------------
% vertical after jump (during shock)
averagingA = [];
averagingM = [];
for (ii=shock1:shock2)
    eA = A(ii,:)/norm(A(ii,:));
    eM = Mnormalized(ii,:)/norm(Mnormalized(ii,:));
    angleAM = acos(eA  *eM')*180/3.14;
    if (abs(angleAM - anglebefore_) < 30)
        averagingA = [averagingA; eA];
        averagingM = [averagingM; eM];
    end;
end;

if (length(averagingA) >= 1)
    if (length(averagingA(:,1)) > 1)
        Vertical2 = mean(averagingA)/norm(mean(averagingA)); % vertical after jump
        M2 = mean(averagingM); M2 = M2/norm(M2);
    else 
        Vertical2 = averagingA/norm(averagingA); % vertical after jump
        M2 = averagingM; M2 = M2/norm(M2);
    end;
        
    meanmagafter = M2 - (M2*Vertical2')*Vertical2;  % horizontal projection of magnetic
    meanmagafter = meanmagafter/norm(meanmagafter); 
    
%  angle of flip is an angle between horizontal projections of magnetic
%  before and after jump
  flip_angle = acos(meanmagafter*meanmagbefore') * 180/3.1416
end;

