function A = randn_cov(n,eigen_mean)

q = randn(n,n);
A = q' * diag(abs(eigen_mean+randn(n,1))) * q;