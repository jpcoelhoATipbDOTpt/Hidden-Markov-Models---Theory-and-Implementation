function q=viterbi_continuous(A,B,c)
% Viterbi algorithm
% q=viterbi_continuous(A,B,c)
% m hidden states, n output states and N observations
% A - mxm (state transitions matrix)
% B - nxm (confusion matrix)
% c - 1xm (priors vector)

[m,N]=size(B);

%% Initialization                     
delta=zeros(N,m);phi=zeros(N,m);
tmp=zeros(m,1);q=zeros(N,1);
t=1;
for k=1:m
    delta(t,k)=c(k)*B(k,t);
    phi(t,k)=0;
end
delta(t,:)=delta(t,:)/sum(delta(t,:));

%% Recursion
for t=2:N,
    for k=1:m,
        for l=1:m,
            tmp(l)=delta(t-1,l)*A(l,k)*B(k,t);
        end
        [delta(t,k),phi(t,k)]=max(tmp);
    end
    delta(t,:)=delta(t,:)/sum(delta(t,:));
end

%% Path finding
[~,Inx]=max(delta(N,:));
q(N)=Inx;
for k=N-1:-1:1,
    q(k)=phi(k+1,q(k+1));
end