function [ACC, NMI, PUR, H, y] = MLPL(KH, HV, truthY, n, m, k, c, alpha, lambda, niter, VN)

%% Initialize parameter
beta = ones(m,1)/sqrt(m);

%% Initialize S and F
S = (KH(:,:,VN)+KH(:,:,VN)')/2;

% initialize F
D = diag(sum(S));
L = D - S;
F = my_eig(L, k, 0);

%% Initialize W
W = zeros(k,k,m);
for v = 1:m
    W(:,:,v) = eye(k);
end

for iter = 1:niter
%     fprintf('Iteration %d...\n', iter);

    %% Upate H
    H = UpdateH(L, W, HV, beta, n, m, k);
    
    %% Update W
    W = UpdateW(HV, H, beta, m, k);
    
    %% Update beta
    beta = Updatebeta(H, HV, W, m);
    
    %% Update HV
    HV = UpdateHV(H, W, beta, m);
    
    %% Update S
    S = UpdateS(H, F, alpha, lambda, c);
    
    
    %% Update F
    S = (S+S')/2;                                                       
    D = diag(sum(S));
    L = D-S;
    F_old = F;
    [F, ~, ev]=my_eig(L, k, 0);
    evs(:,iter+1) = ev;
    
    %update lambda
    thre = 1*10^-10;
    fn1 = sum(ev(1:k));                                             
    fn2 = sum(ev(1:k+1));
    if fn1 > thre
        lambda = 2*lambda;
    elseif fn2 < thre
        lambda = lambda/2;  F = F_old;
    else
        break;
    end
    
end

[clusternum, y] = graphconncomp(sparse(S)); 
y = y';
if clusternum ~= k
    sprintf('Can not find the correct cluster number: %d', clusternum)
end
% y = kmeans(F,k,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
[ACC, NMI, PUR] = ClusteringMeasure(truthY, y);
end

