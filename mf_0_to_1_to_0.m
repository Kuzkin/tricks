%  mf 0_to_1_to_0

	function mf = mf_0_to_1_to_0(X,xmin,xcenter,xmax)
	
		if ((X <= xmin) | (X >= xmax)) 
			mf = 0;
		else
			if (X<=xcenter)
				mf = (X-xmin)/(xcenter-xmin);
			else
				mf = (xmax-X)/(xmax-xcenter);
			end
		end	
				
	