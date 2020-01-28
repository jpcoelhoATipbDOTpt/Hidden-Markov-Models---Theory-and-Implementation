clear all;
clc;

m=2;                 % Number of hidden states
d=3;                 % Observations dimension
Lambda=randn(3,10);  % Observations (taken from a
                     % random process with zero
                     % mean and unity variance)
W=[0.6 0.4;0.5 0.5]; % Weight vector
A=[0.7 0.3;0.2 0.8]; % State transition matrix
c=[0.5;0.5];         % Initial distribution
                     %
mu(1,1,:)=[-1;1;1];  % The Gaussian's parameters
mu(1,2,:)=[2;-2;0];  %
mu(2,1,:)=[-2;-1;0]; %
mu(2,2,:)=[0;-1;0];  %

Sigma(1,1,:,:)=diag([0.2 0.3 0.1]);
Sigma(1,2,:,:)=diag([1 0.5 2]);
Sigma(2,1,:,:)=diag([0.1 0.1 0.1]);
Sigma(2,2,:,:)=diag([1 2 3]);

% Computation of the observations probabilities
B=phi(Lambda,mu,Sigma,W); 

% Forward Algorithm
[Alpha,LogLik] = forward_continuous_norm(A,B,c);

% Backward Algorithm
Beta = backward_continuous_norm(A,B);

% Viterbi Algorithm
q = viterbi_continuous(A,B,c);