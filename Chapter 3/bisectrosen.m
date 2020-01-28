function eta=bisectrosen(gradfunc,x,dk,limI,limS,maxiter)
% Bisection method

for k=1:maxiter,
    eta=.5*(limI+limS);
    x_k=x+eta*dk;
    r=feval(gradfunc,x_k);
    s=r'*dk;
    if s==0,
        break;
    elseif s>0,
        limS=eta;
    else
        limI=eta;
    end
end