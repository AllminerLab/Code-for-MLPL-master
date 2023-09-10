function beta = Updatebeta(H, HV, W, m)
coef = zeros(1,m);
for v=1:m
    coef(1,v) = trace(H'*HV(:,:,v)* W(:,:,v)); 
end    
    beta = coef/norm(coef,2);
end

