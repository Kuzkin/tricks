% find_sitdowns.m
% exclude the case when the set starts from the hole

istart = 1; 
while (S0(istart)>0) %
    istart  = istart+1;
end;

i1 = 1;
i = istart; %i1+1;
jm = 1;
JH = []; % array containing all startsa and ends of holes
 
% Determine JH - array containing all startsa and ends of holes
while (i<N)
	if S0(i)==0 & S0(i+1)>0		% start of local hole
		JH(jm,1) = i+1;   
    end
    if S0(i)>0 & S0(i+1)==0		% end of local hole
		JH(jm,2) = i;
        jm = jm+1;
    end
	i = i+1;
end;

Ssit = zeros(1,N);
ks = 1;

% find middles of sitdowns and their lengths

sitlength = zeros(size(JH)); 
sitmidle  = zeros(size(JH));

for kk = 1:length(JH(:,1))
    sl = JH(kk,2) - JH(kk,1);
    e1 = JH(kk,2);
    e2 = JH(kk,1);
    if (sl>0)    
        Nsmalp = length(find( Anrm(JH(kk,1):JH(kk,2)) < smallacc_smz));
        if (sl > 4 & Nsmalp > 1) % the sitdown is long enough 
            sitmidle(ks)  = round(0.5*(JH(kk,2) + JH(kk,1)));
            sitlength(ks) = sl * h;
            ks = ks + 1;
            for (kkk = JH(kk,1):JH(kk,2))
                Ssit(kkk) = 1;
            end;
        end;
    end;
      
end;

% Ssit is sitdown mf. If it is nonzero then the point belongs to one of sitdowns 
%
% sitmidle is an array containing middle points of all holes. It is used in
% calculate sitdowns.
%
% sitlength is an array containing lengths of holes


ii=1:N;
figure(11); plot(ii, 10*Ssit,'m', ii, Anrm,'b');%, Tg+T0g, V*0.1);


