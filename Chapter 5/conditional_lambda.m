function lambda=conditional_lambda(H,alfa,m)

lambda=zeros(m,1);
for k=1:m,
    lambda(k)=H*alfa(:,k);
end