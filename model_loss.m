function f = model_loss(C, W, params)

% Assume that A is k x n 
%             X is N x k
%               therefore XA is N x n 
%             T is N x n 
%             W is N x n

T = params.optim.T; 
x = params.optim.xpts; 
n = params.optim.n; 
    
R = T*C - repmat(x, 1, n); 
RR = R.*R; 
WR = W.*R; 
GA = T'*WR; 
GW = 0.5*RR; 

f = 0.5*sum(sum(W.*RR));


end

