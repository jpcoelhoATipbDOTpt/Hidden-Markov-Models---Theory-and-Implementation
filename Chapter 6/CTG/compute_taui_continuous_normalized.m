function taui=compute_taui_continuous_normalized(Gamma,B)

[m,N]=size(B);           
taui=Gamma(1:N-1,:);
taui=(sum(taui,1)).'*ones(1,m);