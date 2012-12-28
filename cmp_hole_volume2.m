% cmp_hole_volume2
% compute how much "water" could be held in a hole

function [volume,jmax] = cmp_hole_volume2(Anrm,Tw,hShk,i1,i2,m)

% hole in from i1 to i2
% get boarders as max left and right at within m points from the hole edges

% h1 - hole hight at left
% h2 - hole hight at right (where shock is)

   MAX_LEFT_SEARCH = 10;     % dont search wall farther than that

   [N,mm] = size(Anrm);
   h1 = 0;  h2 = 0; jmax = 0;% jmax - coord of hole wall
   i_left_limit = max([i1-MAX_LEFT_SEARCH,1]);   % dont search left of this point
   % find left edge h1
   hole_max = max(Anrm(i1:i2));     % max Anrm inside the hole
   if i1-m>0 h1 = max(Anrm((i1-m):(i1-1))); end
   % now find local max of left wall
   i = i1-m;
   while i>i_left_limit
     i = i-1;
	 if Anrm(i)>=h1 
	    h1 = Anrm(i);
	 else
		break;
	 end
   end	 
  
   % if left wall top is below hole max - walk left more
   % this may happen when there are trick inside jump but the left end is
   % defined poorely
   if h1<hole_max     % walk left until find larger 
	 h1 = max(Anrm(i_left_limit:(i1-1)));
   end
   jmax = find(Anrm(i_left_limit:(i1-1))==h1);   % for debugging
   % now left wall is the highest that we can find m points left from the
   % hole; 
   % if it is still smaller than max inside the hole - there is no jump
   
     if h1 < hole_max
         volume = 0;
         return;
     end
   
   % for right edge use max of Shock that is passed inside function 
   h2 = hShk;
   h = min([h1,h2]); 
   %s1 = 0.5*(h1+h2)*(i2-i1+1);
   dt = Tw/(i2-i1);
   s1=h*Tw;
   s2 = sum(Anrm(i1:i2))-0.5*(Anrm(i1)+Anrm(i2));
   volume = s1-s2*dt;
  
   