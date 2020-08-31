function LL = loglik(O,HMM)

%% Compute B
B = phi(O,HMM.Mu,HMM.Sigma,HMM.W);
%% Compute log likelihood using the forward algorithm
[~,LL]=forward_continuous_normalized(HMM.A,B,HMM.c);