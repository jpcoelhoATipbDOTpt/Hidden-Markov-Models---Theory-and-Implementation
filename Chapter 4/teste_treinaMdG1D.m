function [sol_W,sol_Mu,sol_Sigma,Fit]=test_trainMdG1D
% Problem 1D @ 2MoG

% (c)2010 João Paulo Coelho

clc
dx=.01;
xx=-1:dx:5;
m=[1 2];
s=[.1 .2];
d=[0.4 0.6];
pdf=dx*(d(1)/(sqrt(2*pi*s(1)))*exp(-0.5*((xx-m(1)).^2)/s(1))+d(2)/(sqrt(2*pi*s(2)))*exp(-0.5*((xx-m(2)).^2)/s(2)));
sum(pdf)
% Generate with this pdf a training vector:
N=1000;
for k=1:N,
    inx=rand;
    if inx<d(1);
        x(k)=m(1)+sqrt(s(1))*randn;
    else
        x(k)=m(2)+sqrt(s(2))*randn;
    end
end
figure;
Nb=100;
[a,b]=hist(x,Nb);
a=(Nb/4)*dx*a/N; % I do not understand why the multiplicative factor...
% show that, in fact, the data have a distribution similar to that defined
plot(b,a,xx,pdf)
% Use of EM to estimate the Gaussian mixture
Niter=1000;
M=2; % Gaussian number
P=zeros(M,N);
p=zeros(1,N);
%
Fit=zeros(Niter,1);
sol_W=zeros(Niter,2);
sol_Mu=zeros(Niter,2);
sol_Sigma=zeros(Niter,2);
% Set initial values
Sigma=3*rand(1,2);
Mu=[-1-2*randn 1+2*randn];
W=rand(1,2);W=W/sum(W);

%-------------------------------------------------------Training Algorithm
[sol_W,sol_Mu,sol_Sigma]=TrainMdG(x,2,100);
%--------------------------------------------------------------------------
plotEM(d,m,s,sol_W,sol_Mu,sol_Sigma)

function plotEM(W,Mu,Sigma,sol_W,sol_Mu,sol_Sigma)

[xx,pdf]=drawGauss(W,Mu,Sigma);
figure
subplot(2,1,1);
[xx,pdf2]=drawGauss(sol_W(1,:),sol_Mu(1,:),sol_Sigma(1,:));
plot(xx,pdf,xx,pdf2);
subplot(2,1,2);
[xx,pdf2]=drawGauss(sol_W(end,:),sol_Mu(end,:),sol_Sigma(end,:));
plot(xx,pdf,xx,pdf2);

function [xx,pdf]=drawGauss(W,Mu,Sigma)
dx=.01;
xx=-1:dx:5;
pdf=(W(1)/(sqrt(2*pi*Sigma(1))))*exp(-0.5*((xx-Mu(1)).^2)/Sigma(1))+(W(2)/(sqrt(2*pi*Sigma(2))))*exp(-0.5*((xx-Mu(2)).^2)/Sigma(2));
