function [KH,alpha] = MyGenerateKerMatrix(X, c)
M = length(X);
N = size(X{1},1);
for v = 1:M
    for i = 1:N
        normItem = std(X{v}(i,:));
        if (0 == normItem)
           normItem = eps;
        end
        X{v}(i,:) = (X{v}(i,:) - mean(X{v}(i,:)))/normItem;
    end
    X{v} = X{v}';
end

for v = 1:M
    distX = L2_distance_1(X{v}, X{v});
    [distX1, idx] = sort(distX, 2);
    KH(:,:,v) = zeros(N);
    rr(:,v) = zeros(N,1);
    for i = 1:N
        di = distX1(i,2:c+2);
        rr(i,v) = 0.5*(c*di(c+1)-sum(di(1:c)));
        id = idx(i,2:c+2);
        KH(i,id,v) = (di(c+1)-di)/(c*di(c+1) - sum(di(1:c)) + eps);
    end    
end
alpha = mean(rr,1);
end

