Magnetic  = Mnormalized(Jump_start2(1):Jump_start2(1)+nL1-1,:); % 

% Tin = 1:0.02:max(Tg);
% GYRin = interp1(Tg, GYR, Tin);
% MAGin = interp1(Tm, MAG, Tin);


for (ig = 1:21)
    M0 = Magnetic(1,:);
    Ggain(ig) = 0.85 + 0.01* (ig - 1);
    L_calib = just_integrate (Ggain(ig)*gyro, [1 0 0 0], h);
   
    dM(ig) =0;
    for (ii=1:nL1)  
        M_estimate(ii,:) = qrot_v3((L_calib(ii,:)), M0);
        dM(ig) =  dM(ig) + norm(M_estimate(ii,:) - Magnetic(ii,:)); 
    end;
end;

figure(5); plot(Ggain, dM)