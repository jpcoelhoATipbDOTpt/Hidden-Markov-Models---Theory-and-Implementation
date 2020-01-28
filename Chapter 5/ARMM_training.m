function [A,c,alfa,sigma2] = ARMM_training(Lambda,m,p,MaxIter)

% Parameters initialization
c = rand(m,1);     % Priors vector
c = c/sum(c);      %
A = rand(m,m);     % Prob. state transition matrix
A = A./(sum(A,2)*ones(1,m));
sigma2 = var(Lambda); % Variance
alfa = zeros(p+1,m);  % Initialize autoreg. param.

for i=1:m,            % Stability is guaranteed
    poles=rand(1,p);                %
    eqx=poly(poles);                %
    alfa(:,i)=[rand -eqx(2:end)].'; %
end                                 %

% Iteration...
for iter=1:MaxIter,
    [P1,P2,~]=...
    hamilton_algorithm(A,c,Lambda,alfa,sigma2);
    P3=kim_algorithm(A,P1,P2,p);
    [alfa,sigma2]=markov_ar_coef(p,Lambda,P3);
    [A,c]=markov_ar_matx(p,A,P3);
end