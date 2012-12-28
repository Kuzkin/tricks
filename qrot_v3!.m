% r = qrot_v3(L,v)

function r = qrot_v3(L,v)
   rr = qrot_v(L,v);
   r = rr(2:4);