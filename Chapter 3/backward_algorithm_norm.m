function Beta=backward_algorithm_norm(A,B,O)
% Backward algorithm with normalization
%[Alfa,LogLik]=algoritmo_directo(A,B,O,c)
% m hidden states, n output states and N observations
% A - mxm (state transitions matrix)
% B - nxm (confusion matrix)
% O - 1xN (observations vector)

[m,~]=size(B);
N=length(O);

%% Initialization
Beta=zeros(N,m);
for k=1:m
	Beta(N,k)=1;
end
v(N)=1/sum(Beta(N,:));      % Scaling coefficient

%% Recursion
for t=N-1:-1:1,
	for i=1:m,
		Beta(t,i)=0;
		for j=1:m,
			Beta(t,i)=Beta(t,i)+A(i,j)*B(j,O(t+1))*Beta(t+1,j);
		end
	end
	v(t)=1/sum(Beta(t,:));  % Scaling coefficient
	Beta(t,:)=v(t)*Beta(t,:);
end