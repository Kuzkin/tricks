% plt_jump_mag.m

% compute range to plot Anrm around jump
% plot data in the range i1:i2 around jump i01:i02
   i1 = max([1,Jump_start2(j)-100]);
   i2 = min([N,Jump_end(j)+100]);
 
% find shock points 
   for (kk = 1:length(Jump_end))
       if (Jump_end(kk) == Jump_end(j))
           jshk = kk; 
           break;
       end;
   end;
   ishk1 = Shk_start(jshk);
   ishk2 = Shk_end(jshk);  
   ttshk1 = Ta(ishk1:ishk2) - Ta(i1);

 figure(44);
   tt1 = Ta(i1:i2) - Ta(i1);
   tt = Ta(ii) - Ta(i1);
   Ta_start = Ta(i1)
   j1 = max([i01-10,1]);
   j2 = min([i02+10,N]);
   ii1 = j1:j2;  % show individual points around jump
   
   if MAG_DATA
      plot(tt1,Anrm(i1:i2), tt1, Mnormalized(i1:i2,:)*5, tt1, Mnrm(i1:i2)*5, tt,Anrm(ii),'.g',tt,avr(ii),'r', ttshk1, Anrm(ishk1:ishk2), '.y')
   else
     plot(tt1,Anrm(i1:i2), Ta(ii1)-Ta(i1),Anrm(ii1),'b.',...   % 
        	 tt,Anrm(ii),'.g',tt,avr(ii),'r', ....             % hole in green, avel hole level - red
			 Ta(ih)-Ta(i1),Anrm(ih),'.r',...                   % left hole edge
			 ttshk1, Anrm(ishk1:ishk2), '.y')                  % shock in yellow
   end
   title(['Jump#',num2str(j),' ',num2str(i01),'-',num2str(i02),' MF=',num2str(mf),...
         ' W=',num2str(ww), ...
         ' Anew/Aout=',num2str(10*r_in_out(k)),...
		 ' Vol=',num2str(fix(vol(k))),...
         ' Vol/W=',num2str(fix(vol(k)/ww))]);
  