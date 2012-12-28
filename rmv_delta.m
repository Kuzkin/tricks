% rmv_delta.m
% remove delta like outliers
% moves array Anrm into array B

Peak_min = 40;
Shldr_max = 30;
Delta_min = 25;

N = length(Anrm);

B = Anrm;

for k=2:(N-1)
   if Anrm(k)>Peak_min & ...
         (Anrm(k)-Anrm(k-1)) > Delta_min & ...
         (Anrm(k)-Anrm(k+1)) > Delta_min & ...
         Anrm(k-1)<Shldr_max & Anrm(k+1)<Shldr_max
      B(k) = (Anrm(k-1)+Anrm(k+1))/2;
   else
      B(k) = Anrm(k);
   end
end
B(N) = Anrm(N);