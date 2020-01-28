function taui=compute_taui(Gama,B,O)
% Compute nu ..
% m hidden states, n output states and N observations
%
% Gama - Nxm (from the forward algorithm)
% O - 1xN (observations vector)
% nu - Return an mxm matrix

[m,~]=size(B);
N=length(O);
taui=Gama(1:N-1,:);
taui=(sum(taui,1)).'*ones(1,m);