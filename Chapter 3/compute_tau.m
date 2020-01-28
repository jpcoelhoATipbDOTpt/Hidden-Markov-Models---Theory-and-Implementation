function tau=compute_tau(Alfa,Beta,A,B,O)
% Compute tau...
% m hidden states, n output states and N observations
% Alfa - Nxm (from the forward algorithm)
% Beta - Nxm (from the backward algorithm)
% A - mxm (state transitions matrix)
% B - nxm (confusion matrix)
% O - 1xN (observations vector)

[m,~]=size(B);
N=length(O);
tau=zeros(m,m);
for k=1:N-1,
    num=A.*(Alfa(k,:).'*Beta(k+1,:)).*(B(:,O(k+1))*ones(1,m)).';
    den=ones(1,m)*num*ones(m,1);
    tau=tau+num/den;
end
