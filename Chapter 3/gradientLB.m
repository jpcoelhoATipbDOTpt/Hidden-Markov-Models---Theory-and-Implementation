function [ALFA,OMEGA]=gradientLB(A,B,O,c)
% Likelihood gradient with respect to B
% m hidden states, n output states and N observations
% A - mxm (state transitions matrix)
% B - nxm (confusion matrix)
% O - 1xN (observations vector)
% c - 1xm (priors vector)
% ALFA - Matrix of direct probabilities
% OMEGA - Matrix of the derivative of L in order to B

[m,n]=size(B);           % n - Number of visible states
N=length(O);             % N - Dimension of the observed sequence
PHI=zeros(m,n,m);
ALFA=zeros(N,m);
OMEGA=zeros(m,n);
PHI2=zeros(m,n,m);
k=1;

for i=1:m,
    ALFA(k,i)=B(i,O(k))*c(i);
    for csi=1:m,
        for delta=1:n,
            if O(k)==delta,
                if i==csi,
                    PHI(csi,delta,i)=c(i);
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
        %----------------------------Computes ALFA----------
        S=0;
        for j=1:m,
            S=S+A(j,i)*ALFA(k-1,j);
        end
        %---------------------------------------------------
        ALFA(k,i)=B(i,O(k))*S;
        for csi=1:m,
            for delta=1:n,
                if O(k)==delta,
                    if i==csi,
                        S1=0;
                        S2=0;
                        for j=1:m,
                            S1=S1+A(j,csi)*ALFA(k-1,j);
                            S2=S2+A(j,csi)*PHI(csi,delta,j);
                        end
                        PHI2(csi,delta,i)=S1+B(csi,O(delta))*S2;
                    else
                        S2=0;
                        for j=1:m,
                            S2=S2+A(j,i)*PHI(csi,delta,j);
                        end
                        PHI2(csi,delta,i)=B(i,O(delta))*S2;
                    end
                else
                    S2=0;
                    for j=1:m,
                        S2=S2+A(j,i)*PHI(csi,delta,j);
                    end
                    PHI2(csi,delta,i)=B(i,O(k))*S2;
                end
                OMEGA(csi,delta)=OMEGA(csi,delta)+PHI2(csi,delta,i);
            end
        end
    end
    PHI=PHI2;
end
