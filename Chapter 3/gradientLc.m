function [c,dLdc]=gradientLc(A,B,O,c)
% Likelihood gradient with respect to c
% m hidden states, n output states and N observations
% A - mxm (state transitions matrix)
% B - nxm (confusion matrix)
% O - 1xN (observations vector)
% c - 1xm (priors vector)

[m,~]=size(B);         % n - Number of visible states
dLdc=zeros(1,m);
k=1;
for i=1:m,
    for xi=1:m,
        if i==xi,
            dLdc(i)=B(xi,O(k));
        else
            dLdc(i)=0;
        end
    end
end