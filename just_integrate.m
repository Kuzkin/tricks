function LRP = just_integrate (G, LZ, h)
  LRP = LZ;
  cur = 1;
  N = length(G);
  
  
  while (cur + 2 <= N)        
     y1 = LRP (cur, :);
     Lf  = Lnext(y1, G, h, cur);
     Lf_in = VectOrth(Lf + LRP(cur,:));
     LRP = [LRP; Lf_in];
     LRP = [LRP; Lf];
     cur = cur + 2;  
  end
  if (cur == N - 1)
    p = length(LRP);
    LRP=[LRP; LRP(p,:)];
  end
end

function Lf = Lnext(y, G, h, num)
   k1 = TransFcn (num, G, y);
   k2 = TransFcn (num + 1, G, y +   h*k1);
   k3 = TransFcn (num + 1, G, y +   h*k2);
   k4 = TransFcn (num + 2, G, y + 2*h*k3);
   Lf  = y  + 2*h*(k1  + 2*k2  + 2*k3  + k4) / 6;
end

function R = TransFcn (n, G, V)
  g1 = G(n, 1);
  g2 = G(n, 2);
  g3 = G(n, 3);
  G = [0; g1; g2; g3];
  R=.5*quatmultiply(G.',V);
end