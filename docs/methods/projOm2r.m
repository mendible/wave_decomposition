function [xout] = projOm2r(xin, k)
% projOm2 projects matrix to nearest matrix with k nnz/row

n = size(xin, 1);
xout = zeros(size(xin));
for i = 1:n
   [~, ind] = sort(abs(xin(i,:)), 'descend');
   xout(i,ind(1:k)) = xin(i, ind(1:k));
end


end

