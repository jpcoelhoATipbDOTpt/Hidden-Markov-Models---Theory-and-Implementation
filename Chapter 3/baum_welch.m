function [A,B,c]=baum_welch(A,B,O,c)
% Baum-Welch algorithm
% m hidden states, n output states and N observations
% A - mxm (state transitions matrix)
% B - nxm (confusion matrix)
% O - 1xN (observations vector)
% c - 1xm (priors vector)
    
Alfa  = forward_algorithm(A,B,O,c);
Beta  = backward_algorithm(A,B,O);
Gama  = compute_gama(Alfa,Beta);
tau   = compute_tau(Alfa,Beta,A,B,O);
taui  = compute_taui(Gama,B,O);
nu    = compute_nu(Gama,B);
Omega = compute_omega(Gama,B,O);
c     = Gama(1,:);
A     = tau./taui;
B     = Omega./nu;