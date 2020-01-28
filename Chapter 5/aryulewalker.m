function alfa = aryulewalker(lambda,p)
mu_lambda = mean(lambda);          % Compute mean
lambda_prime = lambda-mu_lambda;   % Lambda prime
Rox = zeros(1,p+1);                % Initialize rho
for i=0:p,
    Rox(1,i+1) = rho(lambda_prime,0,i);
end
Ro = Rox(2:end).';
P = toeplitz(Rox(1:end-1));        % Toeplitz matrix
alfa = P\Ro;                       % Prm. estimation
alfa = [mu_lambda*(1-sum(alfa));alfa].'; % + dc term 