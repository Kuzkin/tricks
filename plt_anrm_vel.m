% plt_anrm_vel.m

dt = diff(TA);
kk = find(dt>100);
Na = length(TA);
kk=[1;kk;Na];
n = length(kk);
for i=1:(n-1)
   ia1 = kk(i)+1;
   ia2 = kk(i+1)-1;
   iia = ia1:ia2;
   iig = find(T>TA(ia1) & T<TA(ia2));
   plot(TA(iia),Anrm(iia),T(iig),V(iig),'r')
   pause;
end
