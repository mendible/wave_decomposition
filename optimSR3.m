function params = optimSR3(params)

u = params.data.u;
x = params.data.x;
t = params.data.t;

n = params.optim.n;
library = params.optim.library;
mu = params.optim.mu;
zeta = params.optim.zeta;
reg = params.optim.reg; 
lambda = params.optim.lambda;
type = params.optim.type;

%%
dx = x(2)-x(1);
dt = t(2)-t(1);

[Xgrid,Tgrid] = meshgrid(x,t);

%% find the points
[tind,xind] = detect_ridges(u,type);
xpts = x(xind);
tpts = t(tind);

%% setup data
T = library(tpts);
[N,k] = size(T); % number of library functions
X = repmat(xpts,1,n);
params.data.N = N;
pars.T = T;
pars.x = xpts;

loss = @(A,W)lineLoss(A,W,pars);
proxA = @idProx;
proxW = @SimplexProj;

%% k-sparse rows
sparsity = 2;
proxB = @(B,zeta)projOm2c(B, sparsity);
proxBobj = @(x) 0;
lambda_sr3 = 0;

%% l1-prox
%lambda = 0.4;
%proxB = @projL12c;
%proxBobj = @(x)norm(x(:),1);


%% params
sr3Par.reg = reg;
sr3Par.zeta = zeta;
sr3Par.lambda = lambda_sr3;
sr3Par.proxBobj = proxBobj;
sr3Par.mu = mu;
sr3Par.proxB = proxB;
sr3Par.proxW = proxW;
sr3Par.tol = 1e-8;
sr3Par.maxiter = 10000;
sr3Par.X = X;
sr3Par.T = T;
sr3Par.k = k;
sr3Par.n = n;

%% initial separation of poithis i nts
Win = initW(xpts,dx,tpts,dt,n);

%% initial model guess
Ain = initA(xpts,tpts,Win,k,n,lambda,library);
Bin = Ain;
%% run SR3
[Asave, Bsave, Wsave] = SR3(loss, Ain, Bin, Win, sr3Par);

pct = sum(Wsave{end})/N;
L = library(t);
shifts = L(:,1:end-1)*Bsave{end}(1:end-1,:);
shifts = shifts(:,pct>0.01);

params.data.dx = dx;
params.data.dt = dt;

params.SR3.Asave = Asave;
params.SR3.Bsave = Bsave;
params.SR3.Wsave = Wsave;
params.SR3.xpts = xpts;
params.SR3.tpts = tpts;
params.SR3.xind = xind;
params.SR3.tind = tind;
params.SR3.shifts = shifts;

