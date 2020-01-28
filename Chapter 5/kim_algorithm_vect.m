function P3=kim_algorithm_vect(A,P1,P2,p)

[m,N]=size(P1);
P3=zeros(m,N);

% Initialization...
P3(:,N)=P1(:,N);

% Recursion...
for k=N-1:-1:p+1,
    P3(:,k)=P1(:,k).*(A*(P3(:,k+1)./P2(:,k+1)));
end