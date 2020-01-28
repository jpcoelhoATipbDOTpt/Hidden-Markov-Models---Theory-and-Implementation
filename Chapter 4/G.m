function y=G(x,Mu,Sigma)
% Compute uni/multivariate Gaussian
L=size(Sigma);
if length(L)>2,
    Sigma=reshape(Sigma,L(2),L(3));
    y=(1/(sqrt(det(2*pi*Sigma))))*exp(-0.5*(((x-Mu)'/(Sigma))*(x-Mu)));
else
    y=(1/(sqrt(det(2*pi*Sigma))))*exp(-0.5*(((x-Mu)'/(Sigma))*(x-Mu)));
end