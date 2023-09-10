function W = UpdateW(HV, H, beta, m, k)
W = zeros(k,k,m);
for v=1:m
    A = beta(v)*HV(:,:,v)'*H;
    [Uv,~,Vv] = svd(A,'econ');
    W(:,:,v) = Uv*Vv';
end
end

