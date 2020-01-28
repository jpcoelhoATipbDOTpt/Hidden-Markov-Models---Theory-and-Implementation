function roulette=random_func(prob,N)
% N - Number of elements
% prob - Probabilities distribution: vector 1xn with distinct
%        probability distribution
% roulette - Vector 1xN with values within 1 and n according 
%            to distribution prob
if nargin==1,
    N=1;
end
% roulette wheel where each slice is proportional to 
% the probability
L=length(prob);
prob=round(100*prob);
% roulette - vector with 100 elements whose distribution
% depends on P. For example if P=[0.3 0.7] then:
% roulette = [1 1 ... 1 2 2 ...2...2]
%             \-------/ \----------/
%                30          70
roulette=[];
for k=1:L,
    roulette=[roulette k*ones(1,prob(k))];
end

% Generates N values evenly distributed between 1 and 100  
% (it will be the index of the "roulette" vector)
ptr=round(99*rand(1,N)+1);
roulette=roulette(ptr);