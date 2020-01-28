function x=hmmupdate_sol(gradfunc,Theta,O,x,Ai,bi,dk,var)
dk=dk/norm(dk);
% Determine eta max
if ~isempty(Ai)
    den=Ai*dk;
    inx=find(den>0);
    if ~isempty(inx)
        num=bi(inx,:)-Ai(inx,:)*x;
        den=den(inx,:);
        eta=num./den;
        inx=find(eta<0);eta(inx)=[];
        etamax=min(eta);
    else
        etamax=1;
    end
    % It uses the bisection method to determine the optimal for eta:
    eta=hmmbissectrosen(gradfunc,Theta,O,x,dk,0,etamax,100,var);
else
    eta=hmmbissectrosen(gradfunc,Theta,O,x,dk,0,1,100,var);
end
% Calculates new solution
x=x+eta*dk;