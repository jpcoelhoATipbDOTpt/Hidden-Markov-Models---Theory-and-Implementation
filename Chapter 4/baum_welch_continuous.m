function [A,Mu,Sigma,W,c]=baum_welch_cont(O,A,Mu,Sigma,W,c,MaxIter)
% Baum-Welch method for continuous HMM
% [A,Mu,Sigma,W,c]=baum_welch_cont(O,A,Mu,Sigma,W,c,MaxIter)

[D,N]=size(O);[m,M]=size(W);

% Cycle while iterations <= MaxIter
for iter=1:MaxIter
    % Compute B
    B=phi(O,Mu,Sigma,W); %(ok)
    % Compute ALFA using the forward algorithm
    [Alfa,~]=forward_continuous_norm(A,B,c);
    % Iteration info.
    disp(['Iteration - ' num2str(iter) '(Lik = ' num2str(LogLik) ')']);
    % Compute BETA using the backward algorithm
    Beta=backward_continuous_norm(A,B);
    % Compute Gama
    Gamma=compute_gamma(Alfa,Beta);
    % Compute Tau
    tau=tau_continuous(Alfa,Beta,A,B);
    % Compute Taui
    taui=taui_continuous(Gamma,B);
    % Estimation of initial probability c
    c=Gamma(1,:);
    % Estimation of probability transition matrix A
    A=tau./taui;
    % Compute Rho
    rho=rho_continuous(O,Gamma,Mu,Sigma,W);
    % Estimation of weights W
    for i=1:m,
        for j=1:M,
            num=0;
            for k=1:N,
                num=num+rho(k,i,j);
            end
            den=0;
            for k=1:N,
                for l=1:M
                    den=den+rho(k,i,l);
                end
            end
            W(i,j)=num/den;
        end
    end
    % Estimation of mean vector Mu
    for i=1:m,
        for j=1:M,
            num=0;
            den=0;
            for k=1:N,
                num=num+rho(k,i,j)*O(:,k);
                den=den+rho(k,i,j);
            end
            Mu(i,j,:)=num/den;
        end
    end
    % Estimation of covariance Sigma
    for i=1:m,
        for j=1:M,
            num=0;
            den=0;
            for k=1:N,
                M1=reshape(Mu(i,j,:),D,1);
                num=num+rho(k,i,j)*(O(:,k)-M1)*(O(:,k)-M1)';
                den=den+rho(k,i,j);
            end
            Sigma(i,j,:,:)=num/den;
        end
    end
    % Add regularization term to Sigma
    id=eye(D,D);e=zeros(m,M,D,D);
    for i=1:m,
        for j=1:M,
            e(i,j,:,:)=0.01*id;
        end
    end
    Sigma=Sigma+e;
end
