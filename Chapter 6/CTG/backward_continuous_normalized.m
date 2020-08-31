function Beta=backward_continuous_normalized(A,B)

[m,N]=size(B);
%% Initialization
Beta=zeros(N,m);
for k=1:m
    Beta(N,k)=1;
end
u(N)=1/sum(Beta(N,:));
%% Recursion
for t=N-1:-1:1,
    for i=1:m,
        Beta(t,i)=0;
        for j=1:m,
            Beta(t,i)=Beta(t,i)+A(i,j)*B(j,t+1)*Beta(t+1,j);
        end
    end
    u(t)=1/sum(Beta(t,:));
    Beta(t,:)=u(t)*Beta(t,:);
end