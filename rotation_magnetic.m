% determine rotation unsing magnetic and shock
%   Avect = [];
%   for (is = Shk_start(jshk):Shk_end(jshk))
%        Avect(is, :) = A(is,:)/norm(A(is,:)); 
%   end;
%   if  (Shk_end(jshk) - Shk_start(jshk) > 0) 
%       av_A(k,:) = mean(Avect(Shk_start(jshk):Shk_end(jshk),:)) 
%   end;
%   av_A(k,:) =   av_A(k,:) /norm(av_A(k,:));
%  
%    %    % magnetic before and after jump
%    if (exist('MAG'))
%        mag_before(k,:)  = mean(Mnormalized(i01-20:i01-5,:))
%        mag_after(k,:)   = mean(Mnormalized(Shk_end(jshk)+10:Shk_end(jshk)+30,:))
%        angle_magmag(k)  = 180*acos(mag_before(k,:) * mag_after(k,:)')/3.14
%        angle_mag_acc(k) = 180*acos(mean(Mnormalized(Shk_start(jshk):Shk_end(jshk),:)) * av_A(k,:)')/3.14
%    end;
%    % e.o. magnetic before and after jump

% % --- \rottable_mag_selfrot720_20120305
%     i1 = 89;
%     i2 = 196;
%     ishk1 = 196;
%       k = 1;
    
% % % --- salto720 05032012
%     i1 = 100;
%     i2 = 225;
%     ishk1 = 225;
%     k = 1;
% % \
% k=1;
% i1 = 172; i2 = length(Mc);
% ishk1 = 340;
% 
% 
%  for (ii=1:length(MAG)) 
%         Mnormalized(ii,:) = MAG(ii,:)/norm(MAG(ii,:)); 
%         Mnrm(ii) = norm(MAG(ii,:)); 
%         Anrm(ii) = norm(ACC(ii,:)) ; 
%  end;
% 
% 

% angle between vertical and mag

% i1 = i1 - 30;

% ----------------
% % mag_rot_360x3
% k=1;
% i1 = 410;
% ishk1 = 1160;
% i2 = 1200;
% ----------------
% % mag_fast_1
% k=1;
% i1 = 552;
% ishk1 = 842;
% i2 = 1200;

% k=1;
% i1 = 1845;
% ishk1 = 2084;
% i2 = 2000;

% k=1;
% i1 = 2790;
% ishk1 = 3095;
% i2 = 4000;
% ----------------
% % mag_fast_2
% k=1;
% i1 = 590;
% ishk1 = 690;
% i2 = 700;

% % mag_fast_3
% k=1;
i1 = 1180;
ishk1 = 1360;
i2 = 1360;

% Jump_start2(1) = 1180;
% Shk_start(1) = 1360;
% Shk_end(1) =  1360;
% 

angle_mag_vert = 145/180*3.14159; 
sin_angle_mag_vert = sin(angle_mag_vert)

% calculate normal to the plane of rotation
 kk=1;
 dphi = [];
 dphi_1 = [];
 for (is=i1:ishk1-1)
    n_rotplane(kk,:) = cross(Mnormalized(is,:), Mnormalized(is+1,:));
    if (norm(n_rotplane(kk,:)) > 0.001) 
       n_rotplane(kk,:) = n_rotplane(kk,:)/norm(n_rotplane(kk,:));  
       kk = kk + 1;
    end; 
 end;
 n_rotplane_ = mean(n_rotplane);
 n_rotplane_ = n_rotplane_/norm(n_rotplane_);
 % calculate ort of the compoment of Vertical lying in the plane of rotation
 for (is=i1:ishk1-1) 
     Vert_projection(is,:) = Mnormalized(is,:) - n_rotplane_* (n_rotplane_*Mnormalized(is,:)');
%    Vert_projection(is,:) = Vert_projection(is,:)/norm(Vert_projection(is,:));
 end;
 
 mean_Vert = mean(Vert_projection(i1:ishk1-1,:));
 for (is=i1:ishk1-1)    Vert_projection(is,:) = Vert_projection(is,:)  - mean_Vert;   end;
 
 mean_Vert = [0 0 0];
 mean_count = 0;
 for (is=i1:ishk1-2)
    direction  = n_rotplane_ * (cross(Vert_projection(is,:),Vert_projection(is+1,:)))';
    if (direction > 0) mean_Vert = mean_Vert + Vert_projection(is,:);  end; 
    mean_count = mean_count + 1;
 end;
 mean_Vert = mean_Vert/mean_count;
    
 for (is=i1:ishk1-1)    Vert_projection(is,:) = Vert_projection(is,:)  - mean_Vert;   end;
        
    % calculate the angle of salto
    for (is = i1:ishk1-2) 
        direction = sign((Vert_projection(is+1,:)-Vert_projection(is,:)) * (cross(n_rotplane_, Vert_projection(is,:)))');
        dphi(is) = direction * 0.5*(norm(Vert_projection(is,:)) + norm(Vert_projection(is+1,:))) * acos(Vert_projection(is,:)*Vert_projection(is+1,:)'/(norm(Vert_projection(is,:))*norm(Vert_projection(is+1,:))))*180/3.14;
    end;
       
    % calculate the angle of salto
    for (is = i1:ishk1-2) 
        direction = sign((Vert_projection(is+1,:)-Vert_projection(is,:)) * (cross(n_rotplane_, Vert_projection(is,:)))');
        dphi_1(is) = direction * acos(Vert_projection(is,:)*Vert_projection(is+1,:)'/(norm(Vert_projection(is,:))*norm(Vert_projection(is+1,:))))*180/3.14;
    end;
    
    Mstd(k) = std(Mnrm(i1:ishk1-1));
    angle_mag(k) = abs(sum(dphi));
    angle_mag_1(k) = abs(sum(dphi_1));
    
%  ----------------  radius - min/max --------------------------
    ex = Vert_projection(i1,:); ex = ex/norm(ex);
    ey = cross(ex,n_rotplane_); ey = ey/norm(ey);
    for (is = i1:ishk1-1)   vx(is)  = ex*Vert_projection(is,:)';   end;
    for (is = i1:ishk1-1)   vy(is)  = ey*Vert_projection(is,:)';   end;
    
    mean_radius_coord(k) = 0.25*( abs(max(vx(i1:ishk1-1)) - min(vx(i1:ishk1-1))) + abs(max(vy(i1:ishk1-1)) - min(vy(i1:ishk1-1))));
    mean_radius(k) = mean(vnorm(Vert_projection(i1:ishk1-1,:)));

%  ----------------  radius - area --------------------------
    ex = Vert_projection(i1,:); ex = ex/norm(ex);
    ey = cross(ex,n_rotplane_); ey = ey/norm(ey);
    for (is = i1:ishk1-1)   vx(is)  = ex*Vert_projection(is,:)';   end;
    for (is = i1:ishk1-1)   vy(is)  = ey*Vert_projection(is,:)';   end;
    
    Area = 0;
    angle_area = 0;
    mean_vy = mean(vy(i1:ishk1-1));
    for (is = i1:ishk1-2)
        direction  = sign(n_rotplane_ * (cross(Vert_projection(is,:),Vert_projection(is+1,:)))');
       % if (direction > 0) 
%             Area = Area + direction * abs((vx(is+1) - vx(is)) * (0.5 *(vy(is+1) + vy(is)) - mean_vy));
              Area = Area + 0.5  * direction * norm(cross(Vert_projection(is,:), Vert_projection(is+1,:)));
              angle_area = angle_area + direction * (0.5*(norm(Vert_projection(is,:)) + norm(Vert_projection(is+1,:)))/mean_radius(k)) *...
                           abs(acos(Vert_projection(is,:)*Vert_projection(is+1,:)'/(norm(Vert_projection(is,:))*norm(Vert_projection(is+1,:)))));
        %end;
    end;
    
    mean_radius_area(k) = sqrt(abs(2*Area/(angle_area)));
%     angle_area = angle_area;
%     Area = Area;
    
 %   N_Vertical = acos(n_rotplane_* av_A')*180/3.14
    
 % ----------------------   piñtures   ----------------------------------    
    figure(40); plot(i1:i2, Mnrm(i1:i2)/mean(Mnrm(i1:i2)));
   % trajectory of the vertical
    figure(41); plot3(Mnormalized(i1:ishk1-1,1), Mnormalized(i1:ishk1-1,2), Mnormalized(i1:ishk1-1,3));
    hold on
    plot3(1, 1, 1);
    plot3(-1, -1, -1);
    hold off
    % trajectory of the vertical on the plane of rotation
    figure(42); plot3(Vert_projection(i1:ishk1-1,1), Vert_projection(i1:ishk1-1,2), Vert_projection(i1:ishk1-1,3));
    hold on
    plot3([0 Vert_projection(i1,1)], [0 Vert_projection(i1,2)], [0 Vert_projection(i1,3)], 'r');
    plot3(1, 1, 1);
    plot3(-1, -1, -1);
    hold off
    
   
    
    