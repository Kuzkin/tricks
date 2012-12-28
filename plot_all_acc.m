% plot_all_acc

Anrm = vnorm(ACC);
N = length(Ta);
dt = diff(Ta);
jj = find(dt>20);
jj1 = [1,jj',N];
M = length(jj1);

zz = zeros(M-1,1);
AAav = zz;
AAstd = zz;
AAmax = zz;
ii1 = zz; ii2=zz;
figure(10);
for k=1:(M-1)
  i1 = jj1(k)+1; i2 = jj1(k+1)-1;
  ii = i1:i2;
  n = (i2-i1)+1;
  if n < 1000
     plot(Anrm(ii));
	 pause;
  else
	 K = fix(n/500);
	 for m=1:K
	   m1 = i1+(m-1)*500;
	   m2 = m1+500;
	   plot(Anrm(m1:m2));
	   pause
	end
	plot(Anrm(m2+1:i2))
	pause
  end
end	   
  
