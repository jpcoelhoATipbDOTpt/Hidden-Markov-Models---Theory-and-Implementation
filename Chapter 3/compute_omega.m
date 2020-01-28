function Omega=compute_omega(Gama,B,O)
% Compute omega...
% m hidden states, n output states and N observations
% Gama - Nxm matrix
% B - nxm (confusion matrix)
% O - 1xN (observations vector)

[m,n]=size(B);
for j=1:n,
    inx=find(O==j);
    if ~isempty(inx),
        Omega(:,j)=sum(Gama(inx,:),1).';
    else
        Omega(:,j)=0*ones(m,1);
    end
end