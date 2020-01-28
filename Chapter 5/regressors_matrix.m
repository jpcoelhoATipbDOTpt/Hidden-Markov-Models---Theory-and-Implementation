function H=regressors_matrix(Lambda_pred,Lambda,k,p,i)

if i<=p
    H=[1 Lambda_pred(k,i-1:-1:1) Lambda(k:-1:k-(p-i))];
else
    H=[1 Lambda_pred(k,i-1:-1:i-p)];
end
