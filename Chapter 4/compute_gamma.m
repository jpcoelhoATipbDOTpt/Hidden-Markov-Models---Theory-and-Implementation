function Gamma=compute_gamma(Alfa,Beta)
% Compute gamma
% Gamma=compute_gamma(Alfa,Beta)

[~,m]=size(Alfa);
P=diag(Alfa*Beta')*ones(1,m);
Gamma=(Alfa.*Beta)./P;