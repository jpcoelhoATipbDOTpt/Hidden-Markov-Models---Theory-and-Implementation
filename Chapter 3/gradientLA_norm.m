function [ALFA,OMEGA]=gradientLA_norm(A,B,O,c)
% Likelihood gradient with respect to A (normalized version)
% m hidden states, n output states and N observations
% A - mxm (state transitions matrix)
% B - nxm (confusion matrix)
% O - 1xN (observations vector)
% c - 1xm (priors vector)
% ALFA - Matrix of forward probabilities
% OMEGA - Matrix of the derivative of L in order to A

[m,~]=size(B);     % n - Number of visible states
N=length(O);       % N - Dimension of the observed sequence
PHI=zeros(m,m,m);
[ALFA,~,n]=forward_algorithm_norm(A,B,O,c);
OMEGA=zeros(m,m);
PHI2=zeros(m,m,m);
k=1;

for i=1:m,
    ALFA(k,i)=n(k)*B(i,O(k))*c(i);
    for csi=1:m,
        for delta=1:m,
            PHI(csi,delta,i)=0;
            OMEGA(csi,delta)=OMEGA(csi,delta)+PHI(csi,delta,i);
        end
    end
end

for k=2:N,
    for i=1:m,
        S=0;
        for j=1:m,
            S=S+A(j,i)*ALFA(k-1,j);
        end
        ALFA(k,i)=n(k)*B(i,O(k))*S;
        for csi=1:m,
            for delta=1:m,
                if i~=delta,
                    S=0;
                    for j=1:m,
                        S=S+A(j,i)*PHI(csi,delta,j);
                    end
                    PHI2(csi,delta,i)=B(i,O(k))*n(k)*S;
                else
                    S=0;
                    for j=1:m,
                        S=S+A(j,i)*PHI(csi,delta,j);
                    end
                    PHI2(csi,delta,i)=B(i,O(k))*n(k)*(ALFA(k-1,csi)+S);
                end
                OMEGA(csi,delta)=OMEGA(csi,delta)+PHI2(csi,delta,i);
            end
        end
    end
    PHI=PHI2;
end