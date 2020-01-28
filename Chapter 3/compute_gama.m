function Gama=compute_gama(Alfa,Beta)
% Compute gamma
% Alfa - Nxm (from the forward algorithm)
% Beta - Nxm (from the backward algorithm)
% Gama - Return an Nxm matrix with the shape:
%  _                        _
% |  gama_1(1) ... gama_1(m) |
% |   ...             ...    |
% |  gama_N(1) ... gama_N(m) |
%  -                        -

[~,m]=size(Alfa);
P=diag(Alfa*Beta')*ones(1,m);
Gama=(Alfa.*Beta)./P;