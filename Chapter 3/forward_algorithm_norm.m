function [Alfa,LogLik]=forward_algorithm_norm(A,B,O,c)
% Forward algorithm with normalization
%[Alfa,LogLik]=forward_algorithm_norm(A,B,O,c)
% m hidden states, n output states and N observations
% A - mxm (state transitions matrix)
% B - nxm (confusion matrix)
% O - 1xN (observations vector)
% c - 1xm (priors vector)
% 
% Alfa - Partial probability matrix P(Lambda_k ^ q_k).
% LogLik - log likelihood of the observation sequence O


[m,~]=size(B);
N=length(O);
u=zeros(1,N);

%% Initialization 
Alfa=zeros(N,m);
for k=1:m
	Alfa(1,k)=c(k)*B(k,O(1));
end
u(1)=1/sum(Alfa(1,:));       % Scaling coefficient
Alfa(1,:)=u(1)*Alfa(1,:);

%% Recursion
for l=2:N,
	for k=1:m,
		S=0;
		for i=1:m,
			S=S+A(i,k)*Alfa(l-1,i);
		end
		Alfa(l,k)=B(k,O(l))*S;
	end
	u(l)=1/(sum(Alfa(l,:))); % Scaling coefficient
	Alfa(l,:)=u(l)*Alfa(l,:);
end

%% Compute Log Likelihood
LogLik=-sum(log(u));