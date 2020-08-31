function [A,Mu,Sigma,W,c]=baum_welch_continuous_normalized(O,A,Mu,Sigma,W,c,MaxIter)

[D,N]=size(O);[m,M]=size(W);
%% Cycle while iterations <= MaxIter
for iter=1:MaxIter
    %% Compute B
    B=phi(O,Mu,Sigma,W); %(ok)
    %% Compute ALFA using the forward algorithm
    [Alfa,LogLik]=forward_continuous_normalized(A,B,c);
    %% Iteration info.
    disp(['Iteration - ' num2str(iter) '(Lik = ' num2str(LogLik) ')']);
    %% Compute BETA using the backward algorithm
    Beta=backward_continuous_normalized(A,B);
    %% Compute Gama
    Gamma=compute_gamma_normalized(Alfa,Beta);
    %% Compute Tau
    tau=compute_tau_continuous_normalized(Alfa,Beta,A,B);
    %% Compute Taui
    taui=compute_taui_continuous_normalized(Gamma,B);
    %% Estimation of initial probability c
    c=Gamma(1,:);
    %% Estimation of probability transtion matrix A
    A=tau./taui;
    %% Compute Rho
    rho=compute_rho_continuous_normalized(O,Gamma,Mu,Sigma,W);
    %% Estimation of weights W
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
    %% Estimation of mean vector Mu
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
    %% Estimation of covariance Sigma
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
    %% Add regulatization term to Sigma
    id=eye(D,D);e=zeros(m,M,D,D);
    for i=1:m,
        for j=1:M,
            e(i,j,:,:)=0.01*id;
        end
    end
    Sigma=Sigma+e;
end