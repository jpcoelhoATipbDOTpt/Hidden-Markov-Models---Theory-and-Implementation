function x=update_sol(gradF,x,Ai,bi,dk)
% eta value by default...
dk=dk/norm(dk);
% Determine eta max
if ~isempty(Ai)
    eta=(bi-Ai*x)./(Ai*dk);
    etamax=min(eta(eta>0));
    % It uses the bisection method to determine the optimal eta
    % min f(eta)
    % s.t. 0<=eta<=etamax
    eta=bisectrosen(gradF,x,dk,0,etamax(1),100);
else
    eta=bisectrosen(gradF,x,dk,0,1,100);
end
% Calculates new solution
x=x+eta*dk;