% norm of a vector
function v = vnorm(x);
L = length(x);
v = zeros(L,1);
for i=1:L
  z = x(i,:);
  v(i) = sqrt(z*z');
end
