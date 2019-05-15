function [f,GA, GW] = lineLoss(A, W, params)

% Assume that A is k x n 
%             X is N x k
%               therefore XA is N x n 
%             T is N x n 
%             W is N x n

T = params.T; 
x = params.x; 
n = size(A,2); 
    
R = T*A - repmat(x, 1, n); 
RR = R.*R; 
WR = W.*R; 
GA = T'*WR; 
GW = 0.5*RR; 

f = 0.5*sum(sum(W.*RR));


end

