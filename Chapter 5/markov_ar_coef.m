function [alfa,sigma2]=markov_ar_coef(p,Lambda,P3)

% Parameters initialization
Lambda=Lambda(:);  % Lambda is a column vector
[m,N]=size(P3);    % States (m) and observ. (N)
alfa=zeros(p+1,m); % Pre-allocate space
Sigma=zeros(1,m);  % Pre-allocate space
H=zeros(N-p,p+1);  % Pre-allocate space
H(:,1)=1;

% Regressors Matrix
for i=1:p,
    H(:,i+1)=Lambda(p-i+1:end-i);
end
Y=Lambda(p+1:end);

% Weighted least squares parameters estimation
for i=1:m,
    X=P3(i,p+1:end);
    W=diag(X);
    alfa(:,i)=(H'*W*H)\ (H'*W*Lambda(p+1:end));
end

% Variance estimation
for i=1:m
    E=Y-H*alfa(:,i);
    X=P3(i,p+1:end);
    W=diag(X);
    Sigma(i)=E'*W*E;
end

sigma2=sum(Sigma)/(N-p);