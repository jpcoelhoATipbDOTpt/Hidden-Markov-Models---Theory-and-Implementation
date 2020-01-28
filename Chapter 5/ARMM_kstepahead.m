function [L_pred,L_meas] = ARMM_kstepahead(Lambda,A,c,alfa,sigma2,h)

% Parameters initialization
[p,m] = size(alfa);         % Number of states and 
p = p-1;                    % AR order
N=length(Lambda);           % Number of observations
L_pred=zeros(N,h);     % Predicted obs. matrix
L_meas=zeros(N,h);     % Measured obs. matrix

% Execution (1st step)
[~,P2,~]=...
hamilton_algorithm(A,c,L_pred,alfa,sigma2);
for k=p+1:N-h,
    Upsilon_k_k=P2(:,k).';
    % (2nd step)
    for i=1:h,
        H_k_i=...
        regressors_matrix(L_pred,Lambda,k,p,i);
        vec_lambda=conditional_lambda(H_k_i,alfa,m);
        L_pred(k,i)=Upsilon_k_k*(A^h)*vec_lambda;
        L_meas(k,i)=Lambda(k+i-1);
    end
end