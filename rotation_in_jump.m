%  ----------------   Rotation in jump   --------------------------
% The function calculates angles of self-rotation and salto during the
% jump. 
%
% Input parameters are 
% Jump_start2(1)  - jump start
% Shk_start(1) - shock start (end of the jump)
% Shk_end(1)  - shock end;
% Vert_start  - start point of the dataset, where acc is vertical
% Vert_end  - end point of the dataset, where acc is vertical
%
% Output parameters are
% self_rot1, self_rot2 - calculated angles of self rotation
% angle_of_salto - calculated angle of salto


h = 0.02; 

% =================================================
% % % selfrot360 20122011.mat
% Jump_start2(1) = 224;
% Shk_start(1) = 493;
% Shk_end(1) =  496;
% Vert_start = 1;
% Vert_end = 100;
% % 
% % salto360hand 20122011.mat
% Jump_start2(1) = 163;
% Shk_start(1) = 372;
% Shk_end(1) =  376;

% %salto360hand2 20122011.mat
% Jump_start2(1) = 167;
% Shk_start(1) = 393;
% Shk_end(1) =  397;

% % 
% % % salto360mag
Jump_start2(1) = 240;
Shk_start(1) = 406;
Shk_end(1) =  411;
Vert_start = 1;
Vert_end = 100;


% % % % selfrot360mag
% Jump_start2(1) = 133;
% Shk_start(1) = 326;
% Shk_end(1) =  329;

% % % % % % % % % %%% selfrot180salto360 20022012 THE BEST COMPLEX ROT!
% Jump_start2(1) = 148;
% Shk_start(1) = 369;
% Shk_end(1) =  370;

% % works fine with Vertical = [0 0 1] !!!
% % % % % % % %%% salto360selfrot180 20022012
% Jump_start2(1) = 163;
% Shk_start(1) = 448;
% Shk_end(1) =  449;

% %- rot_jump_phone
% Jump_start2(1) = 430;
% Shk_start(1) = 600;
% Shk_end(1) =  600;
% % %- rot_jump_phone
% Jump_start2(1) = 1180;
% Shk_start(1) = 1360;
% Shk_end(1) =  1360;
% 
% % %- rot_jump_phone2
% Jump_start2(1) = 505;
% Shk_start(1) = 631;
% Shk_end(1) =  631;

% % %- rot_phone_bed
% Jump_start2(1) = 468;
% Shk_start(1) = 552;
% Shk_end(1) =  552;
% Vert_start = 1;
% Vert_end = 100;
% Jump_start2(1) = 1210;
% Shk_start(1) = 1280;
% Shk_end(1) =  1280;
% % %- rot_jump_phone3
% Jump_start2(1) = 1509;
% Shk_start(1) = 1700;
% Shk_end(1) =  1700;
% Jump_start2(1) = 2220;
% Shk_start(1) = 2370;
% Shk_end(1) =  2370;
% % %- rot_jump_phone4
% Jump_start2(1) = 1541;
% Shk_start(1) = 1772;
% Shk_end(1) =  1772;
% Jump_start2(1) = 2284;
% Shk_start(1) = 2496;
% % Shk_end(1) =  2496;
% Jump_start2(1) = 3061;
% Shk_start(1) = 3360;
% Shk_end(1) =  3360;
% Jump_start2(1) = 3370;
% Shk_start(1) = 4003;
% Shk_end(1) =  4003;

% % %- rot_phone_salto_360
% Vert_start = 1;
% Vert_end = 100;
% Jump_start2(1) = 1660;
% Shk_start(1) = 1699;
% Shk_end(1) =  1699;

% Jump_start2(1) = 1272;
% Shk_start(1) = 1300;
% Shk_end(1) =  1300;
%
% % % rot_phone_synch
% Jump_start2(1) = 1000;
% Shk_start(1) = 1082;
% Shk_end(1) =  1082;
% Vert_start = 1;
% Vert_end = 100;
%
% % rot_phone_synchr_long
% Jump_start2(1) = 2684;
% Shk_start(1) = 2800;
% Shk_end(1) =  2800;
% Vert_start = 1;
% Vert_end = 100;

% % rot_phone_legsup 
% Jump_start2(1) = 780;
% Shk_start(1) = 840;
% Shk_end(1) =  840;
% Jump_start2(1) = 1300;
% Shk_start(1) = 1450;
% Shk_end(1) =  1450;
% Jump_start2(1) = 1920;
% Shk_start(1) = 2095;
% Shk_end(1) =  2095;
%  Vertical =   mean(A(1600:1800,:)); 

% % % rot_jump_phone_legsup2 
% Jump_start2(1) = 2920;
% Shk_start(1) = 3020;
% Shk_end(1) =  3020;
% Jump_start2(1) = 4000;
% Shk_start(1) = 4194;
% Shk_end(1) =  4194;
% % plusminus180 
% Jump_start2(1) = 1364;
% Shk_start(1) = 1536;
% Shk_end(1) =  1536;
% Vert_start = 2000;
% Vert_end = 2500;
% % % % jump_legsup_norot
% Jump_start2(1) = 3586;
% Shk_start(1) = 3794;
% Shk_end(1) =  3794;

% % % caLIB_GYRO_GAIN2
% Jump_start2(1) = 855;
% Shk_start(1) = 954;
% Shk_end(1) =  954;
% Vert_start = 1;
% Vert_end = 100;

% Jump_start2(1) = 3254;
% Shk_start(1) = 3370;
% Shk_end(1) =  3370;
% Vert_start = 1;
% Vert_end = 100;
% 
% Jump_start2(1) = 1831;
% Shk_start(1) = 1943;
% Shk_end(1) =  1943;
% Vert_start = 1;
% Vert_end = 100;

% =========================================================================

if (exist('ACC')) A = ACC; end;

%    Assume that shock is vertical. Then local coordinates of the vertical
%    vector are  
    Vertical =   mean(A(Vert_start:Vert_end,:)); 
    Vertical = Vertical/norm(Vertical);
    
    nL1 = Shk_end(1) - Jump_start2(1)+ 1; % length of the jump
   
   % calculate orientation of the sportsman at the end of the jump
    Lend =  FindRot_M_only(Vertical, [0 0 1], [1 0 0]); 
    Lend = Lend/norm(Lend);% quaternion at the end of the jump
    gyro =   G(Jump_start2(1):Shk_end(1),:); 
 

%     calculate_gyro_gain; %  Calculate GYRO GAIN using magnetic % currently not used

% calculate quaternions
    L3 = just_integrate (-gyro, Lend, h); % reverse integration
    for (ii=1:nL1) L(ii,:) = L3(nL1-ii+1,:); end; % sort quaternions
% e.o. reverse integration 
% Thus L determines orientation of the sportsman at any moment of time.



% -------------- Calculate the angle of self rotation (selfrot) ------------
if (M)
% calculate normal to the plane of rotation
    kk=1;
    figure(11);
    plot3(1, 1, 1);        plot3(-1, -1, -1);
    hold on;
    n_M = [];
    for (ii=1:length(M)-1)
       if (norm(M(ii,:)) > 0.0001 && norm(M(ii+1,:)) > 0.0001)
           M1 = M(ii,:)/norm(M(ii,:));
           M2 = M(ii+1,:)/norm(M(ii+1,:));
           if (abs(M1*M2') < 0.995) % angle is larger than 5 degrees 
              n_M(kk,:) = cross(M1, M2);
              n_M(kk,:) = n_M(kk,:)/norm(n_M(kk,:));   
              plot3 ([n_M(kk,1) 0], [n_M(kk,2) 0], [n_M(kk,3) 0], 'b'); 
              kk = kk + 1;
           end;
       end;
    end;
    hold off;
    n_M_ = mean(n_M);
    
    if (norm(n_M_) > 0.1) 
        n_M_ = n_M_/norm(n_M_);
        % calculate ort of the compoment of magnetic lying in the plane of rotation
        for (ii=1:length(M)) M_proj(ii,:) = M(ii,:) - n_M_* (n_M_*M(ii,:)');    end;
        % calculate the angle of rotation of magnetic
        for (is = 1:length(M)-1) 
            direction_M(is) = sign((M_proj(is+1,:) - M_proj(is,:)) * (cross(n_M_, M_proj(is,:)))');
            dphi_M(is) = direction_M(is) * 0.5*(norm(M_proj(is,:)) + norm(M_proj(is+1,:)))...
                                         * acos(M_proj(is,:)*M_proj(is+1,:)'/(norm(M_proj(is,:))*norm(M_proj(is+1,:))))...
                                         * 180/3.14;
        end;

        angle_MAG = sum(dphi_M)
        
        % trajectory of magnetic
        figure(77);  
        plot3(1, 1, 1); 
        hold on;        
        plot3(-1, -1, -1);
        for (ii=1:length(M))
             plot3 ([M(ii,1) 0], [M(ii,2) 0], [M(ii,3) 0], 'b'); 
        end;
        hold off
        
    else
        angle_MAG = 0
    end;
end;

    % -------------- Calculate the angle of self rotation (selfrot) ------------
    % Approach 1
    for (ii=1:nL1) 
         om_grav_sign(ii) = gyro(ii, :) * Vertical';
    end;
    selfrot_1 =  abs(sum(om_grav_sign))*180/3.14 * h 

    % Approach 2
    salto = -100000000;
    maxangle1 = 0;     nangle = 100;
    Vplane  = zeros(1, nangle);
    V = zeros(nangle, 3);

    % Consider all vectors X, orthogonal to the Vertical. Calculate the way of these vectors around Vertical.  
    for (kk = 1:nangle)
        angle = 2 * pi/nangle * kk;
        Langle = [cos(0.5*angle) 0 0 sin(0.5*angle)];
        X0 = quatrotate(Langle, [1 0 0]);  
        X = qrot_v3(quatconj(Lend), X0);   X = X/norm(X);
       
        for (ii=1:nL1) 
           V(ii, :) =  cross(gyro(ii, :), X); % velocity of the end of X
           Vplane(ii) = sign(cross(Vertical, X) * V(ii, :)') * norm(V(ii,:) - Vertical * (V(ii,:)*Vertical'));             
        end;
        Vangle(kk) =  sum(Vplane)*180/3.14 * h;
    end;

    selfrot_2 = abs(max(Vangle)) 

    % ----------------------- caclulate angle of salto -----------------------
    % Calculate global coordinates of the Vertical during the jump
    for (ii=1:nL1)     Vert(ii,:) = qrot_v3(L(ii, :), Vertical);    end;
    
    % calculate normal to the plane of rotation
    kk=1;
    for (ii=1:nL1-1)
       n_rotplane(kk,:) = cross(Vert(ii,:), Vert(ii+1,:));
       if (norm(n_rotplane(kk,:)) > 0.00001) 
          n_rotplane(kk,:) = n_rotplane(kk,:)/norm(n_rotplane(kk,:));  
          kk = kk + 1;
       end; 
    end;
    n_rotplane_ = mean(n_rotplane);
    
    % if norm(n_rotplane_) < eps then the rotation was in two different
    % directions
    if (norm(n_rotplane_) > 0.1) 
        n_rotplane_ = n_rotplane_/norm(n_rotplane_);
        % calculate ort of the compoment of Vertical lying in the plane of rotation
        for (ii=1:nL1-1) Vert_projection(ii,:) = Vert(ii,:) - n_rotplane_* (n_rotplane_*Vert(ii,:)');    end;
        % calculate the angle of salto
        for (is = 1:nL1-2) 
            direction(is) = sign((Vert_projection(is+1,:)-Vert_projection(is,:)) * (cross(n_rotplane_, Vert_projection(is,:)))');
            dnorm(is) = norm(Vert_projection(is+1,:)-Vert_projection(is,:));
            dphi(is) = direction(is) * 0.5*(norm(Vert_projection(is,:)) + norm(Vert_projection(is+1,:)))...
                                        * acos(Vert_projection(is,:)*Vert_projection(is+1,:)'...
                                        /(norm(Vert_projection(is,:))*norm(Vert_projection(is+1,:))))*180/3.14;
        end;

        angle_of_salto = sum(dphi)
        
        % % trajectory of the vertical
        figure(3); plot3(Vert(:,1), Vert(:,2), Vert(:,3));
        hold on
        plot3(1, 1, 1);
        plot3(-1, -1, -1);
        hold off
        
        % trajectory of the vertical on the plane of rotation
        figure(4); plot3(Vert_projection(:,1), Vert_projection(:,2), Vert_projection(:,3));
        hold on
        plot3(1, 1, 1);
        plot3(-1, -1, -1);
        hold off
    else
        angle_of_salto = 0
    end;