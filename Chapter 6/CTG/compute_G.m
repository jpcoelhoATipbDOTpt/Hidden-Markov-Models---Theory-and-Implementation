function G=compute_G(O,Mu,Sigma)

[~,N]=size(O);
[m,M,D,~]=size(Sigma);
G = zeros(m,M,N);

for k=1:M
    for i=1:m
        S1=reshape(Sigma(i,k,:,:),D,D);
        M1=reshape(Mu(i,k,:),D,1);
        for j=1:N,
            G(i,k,j)=(1/(sqrt(det(2*pi*S1))))*...
                     exp(-0.5*(((O(:,j)-M1)'/(S1))*(O(:,j)-M1)));
        end
    end
end
