function y=G2(x,Mu,Sigma)
L=size(Sigma);
Sigma=reshape(Sigma,L(2),L(3));
y=(1/(sqrt(det(2*pi*Sigma))))*...
exp(-0.5*(((x-Mu)'/(Sigma))*(x-Mu)));