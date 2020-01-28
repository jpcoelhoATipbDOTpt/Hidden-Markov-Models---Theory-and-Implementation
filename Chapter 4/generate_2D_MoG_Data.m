function x=generate_2D_MoG_Data(u1,v1,u2,v2,w,N)
x=zeros(2,N);       % Training examples matrix
for k=1:N,
    if rand<w(1);
        x(:,k)=u1+sqrt(v1)*randn(2,1);
    else
        x(:,k)=u2+sqrt(v2)*randn(2,1);
    end
end