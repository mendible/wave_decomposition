function [xout] = projOm2c(xin, k)
% projOm2 projects matrix to nearest matrix with k nnz/col
 
n = size(xin, 2); 
xout = zeros(size(xin));
for i = 1:n
   [~, ind] = sort(abs(xin(:,i)), 'descend');
   xout(ind(1:k),i) = xin(ind(1:k),i);
end
 
 
end