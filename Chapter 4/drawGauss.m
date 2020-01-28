function [xx,pdf]=drawGauss(W,Mu,Sigma)
dx=.01;
xx=-1:dx:5;
pdf=(W(1)/(sqrt(2*pi*Sigma(1))))*exp(-0.5*((xx-Mu(1)).^2)/Sigma(1))+(W(2)/(sqrt(2*pi*Sigma(2))))*exp(-0.5*((xx-Mu(2)).^2)/Sigma(2));