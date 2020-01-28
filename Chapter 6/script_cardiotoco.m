% Clear workspace and load data
clear all;
close all;
clc;

load reducedctgdata.mat 
 
% Define model structure
m = 4;             % Number of hidden states
M = 4;             % Number of Gaussian per mixture
MaxIter = 50;      % Number of iterations

% Model Estimation
normalHMM       = ctg_model(Normal,m,M,MaxIter);
suspectHMM      = ctg_model(Suspect,m,M,MaxIter);
pathologicHMM = ctg_model(Pathological,m,M,MaxIter);

% Model Validation
N=length(VNormal);
ResNormal=zeros(N,2);
for i=1:N,
    LL = [loglik(VNormal(:,i),normalHMM), ...
          loglik(VNormal(:,i),suspectHMM), ...
          loglik(VNormal(:,i),pathologicHMM)];
    [~,inx] = max(LL);
    ResNormal(i,1)=1;ResNormal(i,2)=inx;
end

N=length(VSuspect);
ResSuspect=zeros(N,2);
for i=1:N,
    LL = [loglik(VSuspect(:,i),normalHMM), ...
          loglik(VSuspect(:,i),suspectHMM), ...
          loglik(VSuspect(:,i),pathologicHMM)];
    [~,inx] = max(LL);
    ResSuspect(i,1)=2;ResSuspect(i,2)=inx;
end

N=length(VPathologic);
ResPathological=zeros(N,2);
for i=1:N,
    LL = [loglik(VPathological(:,i),normalHMM),...
          loglik(VPathological(:,i),suspectHMM), ...
          loglik(VPathological(:,i),pathologicHMM)];
    [~,inx] = max(LL);
    ResPathological(i,1)=3;ResPathological(i,2)=inx;
end