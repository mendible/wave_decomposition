function [W] = initW(x,dxi,t,dtau,n)

% this function initializes the W matrix for belonging for the SR3
% optimization
 
N = length(x);

% Euclidian distance
E = zeros(N,N);
for jj = 1:N
    for kk = 1:N
        E(jj,kk) = sqrt( (x(jj)-x(kk)).^2 + (t(jj)-t(kk)).^2 );
    end
end

% adjacency matrix
rad = sqrt(dxi^2+dtau^2);
Adj = E < 5 * rad;
Adj = Adj-diag(diag(Adj));

% degree matrix
Deg = diag(sum(Adj,2));

% laplacian matrix
L = Deg-Adj;
[V, ~] = eig(L);

rng(1);
[clust,~] = kmeans(V(:,1),n,'MaxIter',500);

W = zeros(N,n);
for jj = 1:N
    W(jj,clust(jj)) = 1;
end



end

