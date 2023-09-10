function H = UpdateH(L, W, HV, beta, n, m, k)
U = zeros(n,k);
for v=1:m
    U = U + beta(v)*(HV(:,:,v)*W(:,:,v));
end
P = L-U*U';
[H, ~, ~]=my_eig(P, k, 0);
end

