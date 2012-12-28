function vout = qrot_v( q, v )
%  QUATROTATE Rotate a vector by a quaternion.
%   vout = qrot_v(Q,v)
%  v = (1x3) or (1:4) with v(1)=0;
%  vout = (0,x,y,z)

    [n,m] = size(q);
	if m==1 q = q'; end
	if n~=4 & m~=4 error('wrong q dimenton'); end
	
	[n,m] = size(v);
	if m==1 v = v'; end
	if length(v)==3 v = [0,v]; end
	
	qinv = [q(1),-q(2:4)];
	r1 = qmult(q,v);
    vout = qmult(r1,qinv);
