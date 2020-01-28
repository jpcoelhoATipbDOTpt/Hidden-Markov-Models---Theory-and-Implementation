function [Aa,Ai,bi,Nk]=MatGradRest(A,b,C,x)
% Definition of the matrix with inactive constraints
if ~isempty(A)
    Ai=A(find(A*x-b~=0),:);
    bi=b(find(A*x-b~=0),:);
    % Definition of the matrix with active constraints
    Aa=A(find(A*x-b==0),:);
else
    Aa=[];
    Ai=[];
    bi=[];
end
% Construction of the matrix Nk
Nk=[Aa;C];