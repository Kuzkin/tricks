function r = rmmean (A, n1, n2)
  lena = length(A);
  r = A;
  for i = 1:3
      r(1:n1, i) = centrify(A(1:n1, i));
      r(n1:n2, i) = centrify(A(n1:n2, i));
      r(n2:lena, i) = centrify(A(n2:lena, i));
  end
end