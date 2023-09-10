clear
close all

tic
%% load data
%% Mfeat
load('Mfeat.mat');
X = {fac,fou,mor,kar,pix,zer};
dataname = 'Mfeat';
c = 20; 
lambda = 1;
VN = 5; 
fprintf('Dataset: %s, c = %d, lambda = %d ...\n', dataname, c, lambda);
 
%% Parameters
k = length(unique(y));
m = length(X);
n = size(X{1},1);
niter  = 30;

%% Generate multi-kernel matrixes
[KH,alpha] = MyGenerateKerMatrix(X, c);
alpha = mean(alpha);

opt.disp = 0;

HV = zeros(n,k,m);
for v=1:m 
    KH(:,:,v) = (KH(:,:,v)+KH(:,:,v)')/2;
    [Hv, ~] = eigs(KH(:,:,v), k, 'la', opt);
    HV(:,:,v) = Hv;
end

[ACC, NMI, PUR, H, preY] = MLPL(KH, HV, y, n, m, k, c, alpha, lambda, niter, VN);
fprintf('Evalution ...\n');
t = toc;
fprintf('Time is %fs ...\n', t);
fprintf('ACC = %.3f, NMI = %.3f, PUR = %.3f ...\n', ACC, NMI, PUR);
fprintf('------------------------------------------\n');
fprintf('\n');
