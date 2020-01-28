function [ALFA,OMEGA]=gradientLA(A,B,O,c)
% Likelihood gradient with respect to A
% m hidden states, n output states and N observations
% A - mxm (state transitions matrix)
% B - nxm (confusion matrix)
% O - 1xN (observations vector)
% c - 1xm (priors vector)

[m,n]=size(B);           % n - Number of visible states
N=length(O);             % N - Dimension of the observed sequence
PHI=zeros(m,m,m);
ALFA=zeros(N,m);
OMEGA=zeros(m,m);
PHI2=zeros(m,m,m);
k=1;
for i=1:m,
	ALFA(k,i)=B(i,O(k))*c(i);
	for csi=1:m,
		for epsilon=1:m,
			PHI(csi,epsilon,i)=0;
			OMEGA(csi,epsilon)=OMEGA(csi,epsilon)+PHI(csi,epsilon,i);
		end
	end
end
for k=2:N,
	for i=1:m,
		S=0;
		for j=1:m,
			S=S+A(j,i)*ALFA(k-1,j);
		end
		ALFA(k,i)=B(i,O(k))*S;
		for csi=1:m,
			for epsilon=1:m,
				if i~=epsilon,
					S=0;
					for j=1:m,
						S=S+A(j,i)*PHI(csi,epsilon,j);
					end
					PHI2(csi,epsilon,i)=B(i,O(k))*S;
				else
					S=0;
					for j=1:m,
						S=S+A(j,i)*PHI(csi,epsilon,j);
					end
					PHI2(csi,epsilon,i)=B(i,O(k))*(ALFA(k-1,csi)+S);
				end
				OMEGA(csi,epsilon)=OMEGA(csi,epsilon)+PHI2(csi,epsilon,i);
			end
		end
	end
	PHI=PHI2;
end