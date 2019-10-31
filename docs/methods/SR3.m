function [Csave, Bsave, Wsave] = SR3(loss, Cin,  Bin, Win, params)

tol = params.optim.tol;
maxiter = params.optim.maxiter;
conv = 0;

proxB = params.optim.proxB; % e.g. sparse regularizer on the B
proxW = params.optim.proxW;

zeta = params.optim.zeta;
mu = params.optim.mu;

C = Cin;
W = Win;
B = Bin;

T = params.optim.T;
X = params.optim.X;
k = params.optim.k;
n = params.optim.n;

lambda  = params.optim.lambda;
reg = params.optim.reg;

Ik = eye(k);

TT = T.';
eta = 1/zeta;  % lambda now corresponds to zeta

iter = 1;
Csave{1} = Cin;
Wsave{1} = Win;
Bsave{1} = Bin;

while(~conv && iter < maxiter)
    iter = iter + 1;
    
    % update C
    Cnew = zeros(size(C));
    rhs = TT*(W.*X)+eta*B + reg*C;  % CHANGED
    for col = 1:n
        Cnew(:,col) = (TT*(diag(W(:,col))*T)+(eta+reg)*Ik)\rhs(:,col);
    end
    C = Cnew;
    
    % update W
    obj = loss(C,W);
    Wold = W;
    R = 0.5*(T*C-X).^2;
    W = proxW(W-mu*R, 1);
    
    % update B
    Bold = B;
    B = proxB(C, zeta*lambda);
    
    pct = sum(W)/params.data.N;
    Csave{iter} = C(:,pct>0.1);
    Bsave{iter} = B(:,pct>0.1);
    Wsave{iter} = W(:,pct>0.1);
    
    obj = obj + norm(C(:) - B(:))*eta + lambda*params.optim.proxBobj(B);
    
    err = vecnorm(B(:)-Bold(:))/eta + vecnorm(W(:)-Wold(:))/mu;
    %     err = norm(T*A-X);
    conv = err < tol;
    if ~mod(iter,10)
        fprintf('iter: %d, obj: %7.3e, err: %7.16e\n', iter, obj, err);
    end
end

fprintf('iter: %d, obj: %7.3e, err: %7.3e\n', iter, obj, err);

end

