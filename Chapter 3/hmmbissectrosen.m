function eta=hmmbissectrosen(gradfunc,Theta,O,x,dk,limI,limS,maxiter,var)
[m,n]=size(Theta{var});
for k=1:maxiter,
    eta=.5*(limI+limS);
    x_k=x+eta*dk;
    Theta{var}=reshape(x_k,n,m).';
    [~,dLdV]=feval(gradfunc,Theta{1},Theta{2},O,Theta{3});
    row_dLdV=-reshape(dLdV',m*n,1);
    s=row_dLdV'*dk;
    if s==0,
        break;
    elseif s>0,
        limS=eta;
    else
        limI=eta;
    end
end