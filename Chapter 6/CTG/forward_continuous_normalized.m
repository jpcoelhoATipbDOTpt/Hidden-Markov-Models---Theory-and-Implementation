function [Alfa,LogLik]=forward_continuous_normalized(A,B,c)

[m,N]=size(B);
Alfa=zeros(N,m);
u=zeros(1,N);
%% Initialization 
for k=1:m
    Alfa(1,k)=c(k)*B(k,1);
end
u(1)=1/(sum(Alfa(1,:)));   % Scaling coefficient
Alfa(1,:)=u(1)*Alfa(1,:);
%% Recursion
for l=2:N,
    for k=1:m,
        S=0;
        for i=1:m,
            S=S+A(i,k)*Alfa(l-1,i);
        end
        Alfa(l,k)=B(k,l)*S;
    end
    u(l)=1/(sum(Alfa(l,:)));
    Alfa(l,:)=u(l)*Alfa(l,:);
end
%% Compute Log Likelihood
LogLik=-sum(log(u));
