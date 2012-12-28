function y = VectOrth(x)
    numx = numel(x(:,1));
    for j = 1:numx
       y(j,:) = x(j,:) / norm(x(j,:));
    end
end