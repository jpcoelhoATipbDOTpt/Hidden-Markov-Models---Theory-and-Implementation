function N=gaussian(x,mu,sigma)
% Gaussian function
N=(1/sqrt(2*pi*sigma.^2)).*exp(-0.5*((x-mu)./sigma).^2);