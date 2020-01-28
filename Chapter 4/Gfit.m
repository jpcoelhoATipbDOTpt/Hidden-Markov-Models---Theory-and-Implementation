function pdf=Gfit(x,W,Mu,Sigma)
%Compute probability density function
L=size(Sigma);
[D,M]=size(Mu);
if D>1,
    for k=1:length(x);
        pdf(k)=0;
        for j=1:M
            S1=reshape(Sigma(j,:,:),L(2),L(3));
            M1=Mu(:,j);
            pdf(k)=pdf(k)+W(j)/(sqrt(det(2*pi*S1)))*exp(-0.5*(((x(:,k)-M1)'/S1)*(x(:,k)-M1)));
        end
    end
else
    for k=1:length(x);
        pdf(k)=0;
        for j=1:M
            S1=Sigma(j);
            M1=Mu(j);
            pdf(k)=pdf(k)+W(j)/(sqrt(2*pi*S1))*exp(-0.5*(x(:,k)-M1).^2/S1);
        end
    end
end