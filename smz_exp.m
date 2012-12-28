   function Y = smz_exp(X,a)
% smooth matrix X data with exponential smoothing
% smooth forward and back for each column of matrix X
% note, that start and end values have less smoothing

     Y1 = X;
	 Y2 = X;
	 [N,M] = size(X);
	 for k=2:N
	   Y1(k,:) = a*X(k,:)+(1-a)*Y1(k-1,:);
	   i = N-k+1;
	   Y2(i,:) = a*X(i,:)+(1-a)*Y2(i+1,:);
	 end
	 Y = 0.5*(Y1+Y2);
     
   end