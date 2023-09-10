function S = UpdateS(H, F, alpha, lambda, c)
n = size(H,1);
disth = L2_distance_1(H',H');
distf = L2_distance_1(F',F');
[distXs, idx] = sort((disth+lambda*distf)/2, 2);
S = zeros(n);
for i=1:n                                                         
    idxa0 = idx(i,2:c+1);
    distX = distXs(i,idxa0);
    ad = -distX/(2*alpha);
%     ad = -(disth(i,:)+lambda*distf(i,:))/(2*alpha);
    S(i,idxa0) = EProjSimplex_new(ad);
end
end

