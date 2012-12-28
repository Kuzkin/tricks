% plt_jump_mf_w.m

   figure(3);			% plot mf vs W for all jumps
   if pass_flag(k)			% last test to pass
   plot(x,y,'b',w,z,'b',...   % selection region bounds
      W(kk1),Jump_mf(jj(kk1)),'oy',...   % jumps that passed (w,mf) test in yellow
      W(kk2),Jump_mf(jj(kk2)),'oc',...	 % rejected (w,mf) jumps in cyan
      ww, mf, 'og',...                   % passed all tests  - green 
      ww,Final_mf(k),'+g');				 % passed all tests - green (final mf)
   else
   plot(x,y,'b',w,z,'b',...   % selection region bounds
      W(kk1),Jump_mf(jj(kk1)),'oy',...   % jumps that passed (w,mf) test in yellow
      W(kk2),Jump_mf(jj(kk2)),'oc',...	 % rejected (w,mf) jumps in cyan
      ww, mf, 'og',...                   % passed all tests  - green 
      ww,Final_mf(k),'+c');				 % rejected by Final_mf - red
   end
    
   xlabel('Width')
   ylabel('mf')