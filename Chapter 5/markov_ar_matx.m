function [A_new,c]=markov_ar_matx(p,A,P3)

% Parameters initialization
[m,N]=size(P3);   % states (m) and observations(N)
A_new=zeros(m,m); % Pre-allocate space

% Estimation of the new probabilities transition matrix
for i=1:m
    for j=1:m
        numerator=0;
        for k=p+1:N,
            numerator=numerator+A(i,j)*P3(i,k-1);
        end
        denominator=0;
        for k=p+1:N,
            denominator=denominator+P3(i,k-1);
        end        
        A_new(i,j)=numerator/denominator;
    end
end

% New priors estimation
c=P3(:,p+1);