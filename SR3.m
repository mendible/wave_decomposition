function [Asave, Bsave, Wsave] = SR3(loss, Ain,  Bin, Win, params)

tol = params.tol;
maxiter = params.maxiter;
conv = 0;

%proxA = params.proxA;
proxB = params.proxB; % e.g. sparse regularizer on the B
proxW = params.proxW;

zeta = params.zeta;
mu = params.mu;

A = Ain;
W = Win;
B = Bin;

T = params.T;
X = params.X;
k = params.k;
n = params.n;

lambda  = params.lambda;
reg = params.reg;

Ik = eye(k);

TT = T.';
eta = 1/zeta;  % lambda now corresponds to zeta

iter = 1;
Asave{1} = Ain;
Wsave{1} = Win;
Bsave{1} = Bin;

while(~conv && iter < maxiter)
    iter = iter + 1;
    
    % update A
    Anew = zeros(size(A));
    rhs = TT*(W.*X)+eta*B + reg*A;  % CHANGED
    for col = 1:n
        Anew(:,col) = (TT*(diag(W(:,col))*T)+(eta+reg)*Ik)\rhs(:,col);
    end
    %    A = proxA(Anew, -1); % really don't need this
    A = Anew;
    
    % update W
    [obj, ~, GW] = loss(A,W);
    Wold = W;
    R = 0.5*(T*A-X).^2;
    W = proxW(W-mu*R, 1);
    
    % update B
    Bold = B;
    B = proxB(A, zeta*lambda);
    
    Asave{iter} = A;
    Bsave{iter} = B;
    Wsave{iter} = W;
    
    obj = obj + norm(A(:) - B(:))*eta + lambda*params.proxBobj(B);
    
    err = vecnorm(B(:)-Bold(:))/eta + vecnorm(W(:)-Wold(:))/mu;
    %     err = norm(T*A-X);
    conv = err < tol;
    if ~mod(iter,10)
        fprintf('iter: %d, obj: %7.3e, err: %7.16e\n', iter, obj, err);
    end
end

fprintf('iter: %d, obj: %7.3e, err: %7.3e\n', iter, obj, err);

end

