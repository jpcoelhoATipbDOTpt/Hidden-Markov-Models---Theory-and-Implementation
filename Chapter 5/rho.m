function r = rho(lambda,l,i)
r = sum(lambda(1+l:end-i).*lambda(1+i:end-l))/length(lambda);