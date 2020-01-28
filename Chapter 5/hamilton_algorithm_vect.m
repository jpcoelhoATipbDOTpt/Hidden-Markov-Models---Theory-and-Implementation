function [P1,P2,lik]=hamilton_algorithm_vect(A,c,Lambda,alpha,sigma2)

Lambda=Lambda(:); % Just to ensure that Lambda is a column vector...
[p,m]=size(alpha);
p=p-1;

% Observations number.
N=length(Lambda);

% Space allocations
eta=zeros(m,N);P1=zeros(m,N);P2=zeros(m,N);

% For each observation compute eta_k for k=p+1,...,N
for k=p+1:N
    for i=1:m,
        h_k=[1;Lambda(k-1:-1:k-p)];
        eta(i,k)=(1/(sqrt(2*pi*sigma2)))*...
                exp(-((Lambda(k)-h_k'*alpha(:,i)).^2)/(2*sigma2));
    end
end

% Initialize xi
k=p+1;
for i=1:m,
    P1(i,k)=c(i);
end

% Compute P2
for k=p+1:N-1
    P2(:,k)=(P1(:,k).*eta(:,k))/(ones(1,m)*(P1(:,k).*eta(:,k)));
    P1(:,k+1)=A'*P2(:,k);
end
P2(:,N)=(P1(:,N).*eta(:,N))/(ones(1,m)*(P1(:,N).*eta(:,N)));

% Log-Likelihood...
lik=0;
for k=p+1:N,
    lik=lik+log(ones(1,m)*(P1(:,k).*eta(:,k)));
end