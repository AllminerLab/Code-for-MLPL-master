function HV = UpdateHV(H, W, beta, m)
[n, k] = size(H);
HV = zeros(n,k,m);
for v=1:m
    TP = beta(v)*H*W(:,:,v)';
    [Uv,~,Vv] = svd(TP,'econ');
    HV(:,:,v) = Uv*Vv'; 
end
end

