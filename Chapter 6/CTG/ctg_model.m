function HMM = ctg_model(Data,m,M,MaxIter)
%% ----------------------------------------------------------- Markov Model
[d,~]=size(Data);             % Observations dimension
%% ----------------------------------Hidden Markov Model for class = Normal
%  ---------------------------------------------- Parameters initialization
c=rand(m,1);                    % Priors
c=c/sum(c);                     %
A=rand(m,m);                    % Hidden state transition probabilities
A=A./(sum(A,2)*ones(1,m));      %
W=rand(m,M);                    % Mixture weights
W=W./(ones(m,1)*sum(W));        %
mMu=(ones(M,1)*mean(Data,2).'); % Mean initialization
Mu=zeros(m,M,d);                %
for i=1:m,                      %
    Mu(i,:,:)=mMu+randn(M,d);   %
end                             %
mSigma=cov(Data.');             % Covariance initialization
Sigma=zeros(m,M,d,d);           %
for i=1:m,                      %  
    for j=1:M,                  % 
        Sigma(i,j,:,:)=...      %
         mSigma+randn_cov(d,0); %
    end                         %
end                             % 
%  ---------------------------------------------- Train hidden Markov model
[A,Mu,Sigma,W,c] = ...
    baum_welch_continuous_normalized(Data,A,Mu,Sigma,W,c,MaxIter);
%  ---------------------------------------------- Register model parameters
HMM.A = A;
HMM.Mu = Mu;
HMM.Sigma = Sigma;
HMM.W = W;
HMM.c = c;
