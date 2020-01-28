function rho=rho_continuous(O,Gamma,Mu,Sigma,W)
% Compute rho
% rho=rho_continuous(O,Gamma,Mu,Sigma,W)

[N,~]=size(Gamma); 
[m,M]=size(W);
G=compute_G(O,Mu,Sigma);
rho=zeros(N,m,M);

for k=1:N,
    for i=1:m;
        C=0;
        for l=1:M,
            C=C+W(i,l)*G(i,l,k);
        end
        for j=1:M,
            rho(k,i,j)=Gamma(k,i)*W(i,j)*G(i,j,k)/(C+(C==0));
        end
    end
end
