function [ALFA,OMEGA]=gradientLB_norm(A,B,O,c)
% Likelihood gradient with respect to B (normalized version)
% m hidden states, n output states and N observations
% A - mxm (state transitions matrix)
% B - nxm (confusion matrix)
% O - 1xN (observations vector)
% c - 1xm (priors vector)
% ALFA - Matrix of forward probabilities
% OMEGA - Matrix of the derivative of L in order to B

[m,l]=size(B);      % n - Number of visible states
N=length(O);        % N - Dimension of the observed sequence
PHI=zeros(m,l,m);
[ALFA,~,n]=forward_algorithm_norm(A,B,O,c);
OMEGA=zeros(m,l);
PHI2=zeros(m,l,m);
k=1;

for i=1:m,
    for csi=1:m,
        for delta=1:l,
            if O(k)==delta,
                if i==csi,
                    PHI(csi,delta,i)=n(1)*c(i);
                else
                    PHI(csi,delta,i)=0;
                end
            else
                PHI(csi,delta,i)=0;
            end
            OMEGA(csi,delta)=OMEGA(csi,delta)+PHI(csi,delta,i);
        end
    end
end

for k=2:N,
    for i=1:m,
        for csi=1:m,
            for delta=1:l,
                if i==csi,
                    if O(k)==delta,
                        S1=0;
                        S2=0;
                        for j=1:m,
                            S1=S1+A(j,csi)*ALFA(k-1,j);
                            S2=S2+A(j,csi)*PHI(csi,delta,j);
                        end
                         PHI2(csi,delta,i)=n(k)*S1+B(csi,delta)*n(k)*S2;
                    else
                        S2=0;
                        for j=1:m,
                            S2=S2+A(j,csi)*PHI(csi,delta,j);
                        end
                         PHI2(csi,delta,i)=B(csi,O(k))*n(k)*S2;
                    end
                else
                    S2=0;
                    for j=1:m,
                        S2=S2+A(j,i)*PHI(csi,delta,j);
                    end
                    PHI2(csi,delta,i)=B(i,O(k))*n(k)*S2;
                end
                OMEGA(csi,delta)=OMEGA(csi,delta)+PHI2(csi,delta,i);
            end
        end
    end
    PHI=PHI2;
end