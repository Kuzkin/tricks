%  mf 1_to_0

function mf = mf_1_to_0(X,xmin,xmax)

		if (X <= xmin) mf = 1; end
		if (X >= xmax) mf = 0; end
		if (X<xmax & X>xmin) 
			mf = 1 - (X-xmin)/(xmax-xmin);
		end
