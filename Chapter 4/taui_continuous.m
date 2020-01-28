function taui=taui_continuous(Gamma,B)
% Compute taui
% taui=taui_continuous(Gamma,B)

[m,N]=size(B);           
taui=Gamma(1:N-1,:);
taui=(sum(taui,1)).'*ones(1,m);