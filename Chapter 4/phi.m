function B=phi(Lambda,mu,Sigma,W)
% Compute phi
% B=phi(Lambda,mu,Sigma,W)

[m,M,~,~]=size(Sigma);
[d,N]=size(Lambda);
B=zeros(m,N);

% Iteration along each hidden state
for i=1:m
    % Iteration along each observation
    for j=1:N 
        B(i,j)=0;
        % Iteration along each Gaussian
        for k=1:M 
            S=reshape(Sigma(i,k,:,:),d,d);
            V=reshape(mu(i,k,:),d,1);
            B(i,j)=B(i,j) + W(i,k) *(1/sqrt(det(2 * pi * S))) * ...
              exp(-0.5 * ((Lambda(:,j) - V)'/S) *(Lambda(:,j) - V));
        end
    end
end
