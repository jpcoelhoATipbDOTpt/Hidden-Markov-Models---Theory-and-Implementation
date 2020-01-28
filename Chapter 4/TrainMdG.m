function [sol_W,sol_Mu,sol_Sigma]=TrainMdG(O,M,MaxIter)
% O - Vector with observations
% M - Number of Gaussians
[D,N]=size(O); % Number of observations + problem dimension
P=zeros(M,N);
p=zeros(1,N);
% Initialize data containers
LogLik=zeros(MaxIter,1);
sol_W=zeros(MaxIter,M);
sol_Mu=zeros(MaxIter,D,M);
sol_Sigma=zeros(MaxIter,M,D,D);
Sigma=zeros(M,D,D);
Mu=zeros(M,D);
% ........................................set initial values
for j=1:M
    Sigma(j,:,:)=3*rand(D,D).*eye(D,D);
    Mu(:,j)=randn(D,1);
end
W=rand(1,M);W=W/sum(W);
% .......................................iterative procedure
for iter=1:MaxIter
    LogLik(iter)=sum(log(Gfit(O,W,Mu,Sigma)));
    disp('................................................');
    disp(['Iteration' num2str(iter-1) ': logLik = ' ...
 	      num2str(LogLik(iter))]);
    sol_W(iter,:)=W;
    sol_Mu(iter,:,:)=Mu;
    sol_Sigma(iter,:,:,:)=Sigma;   
    % determine P(Gj|xi)
    for i=1:N,
        p(i)=0;
        for j=1:M
            p(i)=p(i)+W(j)*G(O(:,i),Mu(:,j),Sigma(j,:,:));
        end
        for j=1:M,
            P(j,i)=W(j)*G(O(:,i),Mu(:,j),Sigma(j,:,:))/p(i);
        end
    end
    % compute new estimates for Mu, Sigma and W
    %.....................................................W
    for l=1:M
        soma=0;
        for i=1:N,
            soma=soma+P(l,i);
        end
        W(l)=(1/N)*soma;
    end
    %.....................................................Mu
    for l=1:M
        num=0;
        den=0;
        for i=1:N,
            num=num+O(:,i)*P(l,i);
            den=den+P(l,i);
        end
        
        Mu(:,l)=num./den;
    end
    %...................................................Sigma
    for l=1:M
        num=0;
        den=0;
        for i=1:N,
            num=num+P(l,i)*((O(:,i)-Mu(:,l))*(O(:,i)-Mu(:,l))');
            den=den+P(l,i);
        end
        Sigma(l,:,:)=num./den;
    end
end