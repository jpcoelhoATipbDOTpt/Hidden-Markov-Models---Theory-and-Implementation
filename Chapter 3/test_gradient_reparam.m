function  likelihood=test_gradient_reparam(O,m,n,eta,MaxIter)
% model characteristics...
% m- number of hidden states
% n- number of observable states
% eta - learning rate (vector)
% MaxIter - maximum number of iterations

taA=eta(1); etaB=eta(2); etac=eta(3); % learning rate
% ..........................PARAMETERS INITIALIZATION
W=randn(m,m);A=exp(W)./((sum(exp(W).'))'*ones(1,m));
Q=randn(m,n);B=exp(Q)./((sum(exp(Q).'))'*ones(1,n));
Y=randn(m,1);c=exp(Y)./(sum(exp(Y)));
likelihood=zeros(1,MaxIter);

for j=1:MaxIter
    [~,OMEGA]=gradientLA_norm(A,B,O,c);
    dAdW=A.*(1-A);
    W=W+etaA*OMEGA.*dAdW;
    [~,OMEGA]=gradientLB_norm(A,B,O,c);
    dBdW=B.*(1-B);
    Q=Q+etaB*OMEGA.*dBdW;
    [dLdc]=gradientLc(B,O);
    dcdW=c.*(1-c);
    Y=Y+etac*dLdc.*dcdW;
    A=exp(W)./((sum(exp(W).'))'*ones(1,m));
    B=exp(Q)./((sum(exp(Q).'))'*ones(1,n));
    c=exp(Y)./(sum(exp(Y)));
    [~,lik]=forward_algorithm_norm(A,B,O,c);
    disp(['Iteration -- ' num2str(j) ' FIT -- ' num2str(lik)]);
    likelihood(j)=lik;
end
