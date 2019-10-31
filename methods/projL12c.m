function [xout] = projL12c(xin, step)
% projOm2 projects matrix to nearest matrix with k nnz/col
 
n = size(xin, 2); 
xout = zeros(size(xin));
for i = 1:n
   xout(:,i)  = sign(xin(:,i)).*max(abs(xin(:,i))-step,0); 
end
 
 
end