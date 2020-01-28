function varargout=forward_algorithm(A,B,O,I)
% Forward Algorithm for discrete hidden Markov Models with 'm' hidden
% states, 'n' observable states and 'N' observations.
%   A - mxm (state transition matrix)
%   B - mxn (confusion matrix)
%   O - 1xN (observations vector)
%   c - 1xm (initial probabilities vector)
%
% Usage:
%   [Alfa]=forward_algorithm(A,B,O,I)
%   [Alfa,LogLik]=forward_algorithm(A,B,O,I)
%
% Where:
%   Alfa - Partial probability matrix P(Lambda_k ^ q_k).
%   LogLik - Log-likelihood of O.

[m,n]=size(B);
N=length(O);
 
%% Initialization 
Alfa=zeros(N,m);
for k=1:m
    Alfa(1,k)=I(k)*B(k,O(1));
end

%% Recursion
for l=2:N,
    for k=1:m,
        S=0;
        for i=1:m,
            S=S+A(i,k)*Alfa(l-1,i);
        end
        Alfa(l,k)=B(k,O(l))*S;
    end
end

% Probability of observing O
P=sum(Alfa(N,:));

% function return
if nargout==1,
    varargout={Alfa};
else
    varargout(1)={Alfa};varargout(2)={log(P)};
end