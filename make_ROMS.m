function ROMS = make_ROMS(optim_params,rank,lambda)

u = optim_params.u;
shifts = optim_params.shifts;
t = optim_params.t;
dx = optim_params.dx;

% do POD on unshifted data and truncate
[U.pod,S.pod,V.pod] = svd(u);
upod = U.pod(:,1:rank)*S.pod(1:rank,1:rank)*V.pod(:,1:rank)';

for j = 1:size(shifts,2)
    
    % shift for each wave speed
    for jj = 1:length(t)
        ushift(:,jj,j) = circshift(u(:,jj),-round(shifts(jj,j)/dx));
    end
    
    % do POD on the shifted data
    [U.spod{j},S.spod{j},V.spod{j}] = svd(ushift(:,:,j));
    lowspod(:,:,j) = U.spod{j}(:,1:2*rank)*S.spod{j}(1:2*rank,1:2*rank)*V.spod{j}(:,1:2*rank)';
    
    % do RPCA on the shifted data
    [LowRank{j},Sparse{j},rpcaiter] = inexact_alm_rpca(ushift(:,:,j),lambda);
    [U.srpca{j},S.srpca{j},V.srpca{j}] = svd(LowRank{j});
    lowsrpca(:,:,j) = U.srpca{j}(:,1:rank)*S.srpca{j}(1:rank,1:rank)*V.srpca{j}(:,1:rank)';
    
    % unshift for each wave speed
    for jj = 1:length(t)
        usrpca(:,jj,j) = circshift(lowsrpca(:,jj,j),round(shifts(jj,j)/dx));
        uspod(:,jj,j) = circshift(lowspod(:,jj,j),round(shifts(jj,j)/dx));
    end
    
end


ROMS.U = U;
ROMS.S = S;
ROMS.V = V;
ROMS.usrpca = sum(usrpca,3)/2;
ROMS.uspod = sum(uspod,3)/2;
ROMS.upod = upod;
end
