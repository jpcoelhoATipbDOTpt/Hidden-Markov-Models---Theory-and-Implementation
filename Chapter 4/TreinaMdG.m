function [sol_W,sol_Mu,sol_Sigma]=TrainMdG(O,M,MaxIter)

% O - Vector with observations
% M - Number of Gaussians

% (c)2010 João Paulo Coelho

% Use of EM to estimate the Gaussian mixture

[D,N]=size(O); % Number of observations + Dimension of the problem
P=zeros(M,N);
p=zeros(1,N);

Fit=zeros(MaxIter,1);
sol_W=zeros(MaxIter,M);
sol_Mu=zeros(MaxIter,D,M);
sol_Sigma=zeros(MaxIter,M,D,D);
% ................................................set initial values
for j=1:M
    Sigma(j,:,:)=3*rand(D,D).*eye(D,D);
    Mu(:,j)=[randn(D,1)];
end
W=rand(1,M);W=W/sum(W);

for iter=1:MaxIter
    Fit(iter)=sum(log(Gfit(O,W,Mu,Sigma)));
    disp('..............................................................');
    disp(['Iteration' num2str(iter-1) ': logLik = ' num2str(Fit(iter))]);
    
    sol_W(iter,:)=W;
    sol_Mu(iter,:,:)=Mu;
    sol_Sigma(iter,:,:,:)=Sigma;
    
    % determine P(Gj|xi) - membership probability
    for i=1:N,
        p(i)=0;
        for j=1:M
            p(i)=p(i)+W(j)*G(O(:,i),Mu(:,j),Sigma(j,:,:));
        end
        for j=1:M,
            P(j,i)=W(j)*G(O(:,i),Mu(:,j),Sigma(j,:,:))/p(i);
        end
    end
    % determine new estimates for Mu, Sigma and W
    %............................................................W
    for l=1:M
        soma=0;
        for i=1:N,
            soma=soma+P(l,i);
        end
        W(l)=(1/N)*soma;
    end
    %............................................................Mu
    for l=1:M
        num=0;
        den=0;
        for i=1:N,
            num=num+O(:,i)*P(l,i);
            den=den+P(l,i);
        end
        
        Mu(:,l)=num./den;
    end
    %............................................................Sigma
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
%----------------------------------------------------------------sub-function
function y=G(x,Mu,Sigma)
%(c)2010 João Paulo Coelho
L=size(Sigma);
if length(L)>2,
    Sigma=reshape(Sigma,L(2),L(3));
    y=(1/(sqrt(det(2*pi*Sigma))))*exp(-0.5*(((x-Mu)'/(Sigma))*(x-Mu)));
else
    y=(1/(sqrt(det(2*pi*Sigma))))*exp(-0.5*(((x-Mu)'/(Sigma))*(x-Mu)));
end


%----------------------------------------------------------------sub-function
function pdf=Gfit(x,W,Mu,Sigma)
%(c)2010 João Paulo Coelho
L=size(Sigma);
[D,M]=size(Mu);
if D>1,
    for k=1:length(x);
        pdf(k)=0;
        for j=1:M
            S1=reshape(Sigma(j,:,:),L(2),L(3));
            M1=Mu(:,j);
            pdf(k)=pdf(k)+W(j)/(sqrt(det(2*pi*S1)))*exp(-0.5*(((x(:,k)-M1)'/S1)*(x(:,k)-M1)));
        end
    end
else
    for k=1:length(x);
        pdf(k)=0;
        for j=1:M
            S1=Sigma(j);
            M1=Mu(j);
            pdf(k)=pdf(k)+W(j)/(sqrt(2*pi*S1))*exp(-0.5*(x(:,k)-M1).^2/S1);
        end
    end
end
