% plt_all_jump_mf_w.m

figure(3);			% plot mf vs W for all jumps

kf1 = find(pass_flag == TRUE); % jumps that passed (w,mf) test for Final_mf 
kf2 = find(pass_flag == FALSE); % jumps that have not passed (w,mf) test for Final_mf 

if n_jumps_good > 0
plot(x,y,'b',w,z,'b',...   % selection region bounds
  W(kk1),Jump_mf(i_jump),'or',...   % jumps that passed (w,mf) in red - default they did not pass
  W(kk2),Jump_mf(jj(kk2)),'oc',...	 % rejected (w,mf) jumps in cyan
  ALL_JUMPS(:,3), Jump_mf(k_jumps_good),'og',... % passed all tests - green circle
  widht(kf1), Final_mf(kf1), '+g');		%  Final_mf passed all tests- green cross
else
plot(x,y,'b',w,z,'b',...   % selection region bounds
  W(kk1),Jump_mf(i_jump),'or',...   % jumps that passed (w,mf) in red - default
  W(kk2),Jump_mf(jj(kk2)),'oc',...  	 % rejected (w,mf) jumps in cyan
  widht(kf2), Final_mf(kf2), '+c');		 % Final_mf is rejected - cyan cross
end

xlabel('Width')
ylabel('mf')