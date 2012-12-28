%  mf 0_to_1

	function mf = mf_0_to_1(X,xmin,xmax)
	
		if (X <= xmin) mf = 0; end
		if (X >= xmax) mf = 1; end
		if (X<xmax & X>xmin) 
			mf = (X-xmin)/(xmax-xmin);
		end
			