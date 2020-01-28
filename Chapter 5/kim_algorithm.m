function P3 = kim_algorithm(A,P1,P2,p)

% Parameters initialization
[m,N]=size(P1);
P3=zeros(m,N);
for i=1:m,
    P3(i,N)=P2(i,N);
end

% Recursion
for k=N-1:-1:p+1
    for i=1:m,
        second_term=0;
        for j=1:m
            second_term=second_term+A(i,j)*P3(j,k+1)/P1(j,k+1);
        end
        P3(i,k)=P2(i,k)*second_term;
    end
end