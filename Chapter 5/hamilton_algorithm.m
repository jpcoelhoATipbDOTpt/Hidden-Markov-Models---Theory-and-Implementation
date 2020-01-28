function [P1,P2,lik]=hamilton_algorithm(A,c,Lambda,alfa,sigma2)
% Parameters initialization
Lambda = Lambda(:);  % Lambda is a column vector
[n,m] = size(alfa);  % order (p) and states (m)
p = n-1;             % p = nbr coeff. - 1
N = length(Lambda);  % Number of observations.
lik = 0;             % Model likelihood
P0 = zeros(m,N);  % P(\lambda_k|q_k ^ \Lambda_{k-1})
P1 = zeros(m,N);  % P(q_k=s_i|\Lambda_{k-1})
P2 = zeros(m,N);  % P(q_k=s_i|\Lambda_k)
% Impose P(q_k=s_i|\Lambda_{k-1})=c_i for k = p+1
k = p + 1;
for i=1:m,
    P1(i,k)=c(i);
end
for k=p+1:N-1,
    % Compute P(q_k=s_i|\Lambda_{k}) and 
    % P(\lambda_k|q_k \wedge \Lambda_{k-1})
    for i=1:m,
        h=[1;Lambda(k-1:-1:k-p)];
        P0(i,k)=(1/(sqrt(2*pi*sigma2)))*...
     exp(-((Lambda(k)-h'*alfa(:,i)).^2)/(2*sigma2));
        P2(i,k)=P0(i,k)*P1(i,k);
    end
    denominator=sum(P2(:,k));
    lik=lik+log(denominator);
    for i=1:m,
        P2(i,k)=P2(i,k)/denominator;
    end
    % Compute P(q_k=s_i|\Lambda_{k-1})
    for i=1:m,
        P1(i,k+1)=0;
        for j=1:m,
            P1(i,k+1)=P1(i,k+1)+A(j,i)*P2(j,k);
        end
    end
end
% Compute  P(q_k=s_i|\Lambda_{k}) for k = N
k=N;
for i=1:m,
    h=[1;Lambda(k-1:-1:k-p)];
    P0(i,k)=(1/(sqrt(2*pi*sigma2)))*...
    exp(-((Lambda(k)-h'*alfa(:,i)).^2)/(2*sigma2));
    P2(i,k)=P0(i,k)*P1(i,k);
end
denominator=sum(P2(:,k));
lik=lik+log(denominator);
for i=1:m,
    P2(i,k)=P2(i,k)/denominator;
end