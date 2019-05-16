function params = optimSR3(params)

u = params.data.u;
x = params.data.x;
t = params.data.t;

n = params.optim.n;
library = params.optim.library;
init_lambda = params.optim.init_lambda;
type = params.optim.type;

%%
dx = x(2)-x(1);
dt = t(2)-t(1);

[Xgrid,Tgrid] = meshgrid(x,t);

%% find the points
[tind,xind] = detect_ridges(u,type);
xpts = x(xind);
tpts = t(tind);
params.optim.xind = xind;
params.optim.tind = tind;
params.optim.xpts = xpts;
params.optim.tpts = tpts;

%% setup data
T = library(tpts);
params.optim.T = T;

[N,k] = size(T); % number of library functions
params.data.N = N;
params.optim.k = k;

X = repmat(xpts,1,n);
params.optim.X = X;

loss = @(C,W)model_loss(C,W,params);
proxW = @SimplexProj;
params.optim.proxW = proxW;

%% k-sparse rows
sparsity = 2;
proxB = @(B,zeta)projOm2c(B, sparsity);
params.optim.proxB = proxB;

proxBobj = @(x) 0;
params.optim.proxBobj = proxBobj;

%% l1-prox
%lambda = 0.4;
%proxB = @projL12c;
%proxBobj = @(x)norm(x(:),1);

%% initial separation of poithis i nts
Win = initialize_W(xpts,dx,tpts,dt,n);

%% initial model guess
Cin = initialize_C(xpts,tpts,Win,k,n,init_lambda,library);
Bin = Cin;
%% run SR3
params.optim.lambda = 0;
params.optim.tol = 1e-8;
params.optim.maxiter = 10000;
params.optim.mu = 1e8;

[Csave, Bsave, Wsave] = SR3(loss, Cin, Bin, Win, params);

pct = sum(Wsave{end})/N;
L = library(t);
shifts = L(:,1:end-1)*Bsave{end}(1:end-1,:);
shifts = shifts(:,pct>0.01);

params.data.dx = dx;
params.data.dt = dt;

params.optim.Csave = Csave;
params.optim.Bsave = Bsave;
params.optim.Wsave = Wsave;
params.optim.shifts = shifts;

