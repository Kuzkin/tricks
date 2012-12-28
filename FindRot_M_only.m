function [L] = FindRot_M_only(M, M_init, G_init)
%   The function determines orientation using either magnetic vector only (if use_twitch == 0) or 
%   both magnetic vector and acceleration (at pooints where A~Gravity)
    flag = 0;     
    e = cross(M, M_init);
    if (norm(e) < 0.01) flag = 1; end;
          
    if (flag == 0) 
          e = e/norm(e);      
          kk = e;
          ii = M/norm(M);
          jj = cross(kk, ii);
          jj = jj/norm(jj);

          X = M_init*ii';
          Y = M_init*jj'; 

          phi_temp = atan2(Y,X);  % either X or may be zero      
          if (phi_temp >= 0.001)
              phi_e =  -0.5 * phi_temp;
          else
              phi_e =  -0.5 * (2*pi + phi_temp);
          end;
           e = e * sin(phi_e);
           L_e = [cos(phi_e) e(1,1) e(1,2) e(1,3)];
           L_e = L_e/norm(L_e);
    else 
           if (M_init * M' > 0) 
              L_e = [1 0 0 0];
           else
              m0 = G_init - (M_init* G_init')*M_init/norm(M_init)^2; 
              m0 = m0/norm(m0);
              L_e = [0 m0(1) m0(2) m0(3)];
           end;
    end;                  
    L = L_e;
end

