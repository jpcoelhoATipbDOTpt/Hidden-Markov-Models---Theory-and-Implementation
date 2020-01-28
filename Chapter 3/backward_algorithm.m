function Beta=backward_algorithm(A,B,O)

% Backward Algorithm for discrete hidden Markov Models with 'm' hidden
% states, 'n' observable states and 'N' observations.
%   A - mxm (state transition matrix)
%   B - mxn (confusion matrix)
%   O - 1xN (observations vector)

[m,n]=size(B);
N=length(O);

%% Initialization
Beta=zeros(N,m);
for k=1:m
    Beta(N,k)=1;
end

%% Recursion
for t=N-1:-1:1,
    for i=1:m,
        Beta(t,i)=0;
        for j=1:m,
            Beta(t,i)=Beta(t,i)+A(i,j)*B(j,O(t+1))*Beta(t+1,j);
        end
    end
end