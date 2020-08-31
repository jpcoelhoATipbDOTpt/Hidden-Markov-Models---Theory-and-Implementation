function Gamma=compute_gamma_normalized(Alfa,Beta)

[~,m]=size(Alfa);
P=diag(Alfa*Beta')*ones(1,m);
Gamma=(Alfa.*Beta)./P;