function B=phi(Lambda,mu,Sigma,W)

[m,M,~,~]=size(Sigma);
[d,N]=size(Lambda);
B=zeros(m,N);
for i=1:m %.............................. Iteration along each hidden state
    for j=1:N %........................... Iteration along each observation
        B(i,j)=0;
        for k=1:M %.......................... Iteration along each Gaussian
            S=reshape(Sigma(i,k,:,:),d,d);
            V=reshape(mu(i,k,:),d,1);
            B(i,j)=B(i,j) + W(i,k) * (1/sqrt(det(2 * pi * S))) * ...
                   exp(-0.5 * ((Lambda(:,j) - V)'/S) * (Lambda(:,j) - V));
        end
    end
end