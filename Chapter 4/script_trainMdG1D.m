clear all;
close all;
clc;
%------------------------------------Generate Training Data
dx=.01;
xx=-1:dx:5;
m=[1 2];     % Gaussian Means
s=[.1 .2];   % Gaussian Variances
d=[0.4 0.6]; % Weighting factors
% Generate a training vector with this PDF
N=2000;      % Number of samples
x=zeros(1,N);% Traing samples vector
for k=1:N,
    inx=rand;
    if inx<d(1);
        x(k)=m(1)+sqrt(s(1))*randn;
    else
        x(k)=m(2)+sqrt(s(2))*randn;
    end
end
% Use of EM to estimate the Gaussian mixture
Maxiter=1000;   % Number of iterations
M=2;          % Number of Gaussian functions in the mixture
%----------------------------------------Training Algorithm
[sol_W,sol_Mu,sol_Sigma]=TrainMdG(x,M,Maxiter);
%----------------------------------------------------------
plotEM(d,m,s,sol_W,sol_Mu,sol_Sigma);
