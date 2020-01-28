function [xx,pdf]=drawGauss2(W,Mu,Sigma)
L=size(Sigma);
if length(L)>3
    v1=reshape(Sigma(1,1,:,:),L(3),L(4));
    v2=reshape(Sigma(1,2,:,:),L(3),L(4));
else
    v1=reshape(Sigma(1,:,:),L(2),L(3));
    v2=reshape(Sigma(2,:,:),L(2),L(3));
end
L=size(Mu);
if length(L)>2
    Mu=reshape(Mu(1,:,:),L(2),L(3));
end
dx=.1;
[x1,x2]=meshgrid(-1:dx:4);
xx=[x1(:) x2(:)]';
for k=1:length(xx);
    pdf(k)=W(1)/(sqrt(det(2*pi*v1)))*exp(-0.5*(((xx(:,k)-Mu(:,1))'/v1)*(xx(:,k)-Mu(:,1))))+...
        W(2)/(sqrt(det(2*pi*v2)))*exp(-0.5*(((xx(:,k)-Mu(:,2))'/v2)*(xx(:,k)-Mu(:,2))));
end