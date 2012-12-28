function y = VectNorm (x)
    N = length(x);
    for i = 1:N
        y(i) = norm(x(i,:));
    end
end