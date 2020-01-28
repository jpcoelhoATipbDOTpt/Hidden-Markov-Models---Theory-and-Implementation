% Generate the time-series
clear all;
clc;
y(1)=0;y(2)= 0.06;y(3)=0.2;y(4)=0.4;y(5)=0.5;
for k=5:200,
    y(k)=2.7*y(k-1)-3.1*y(k-2)+1.8*y(k-3)-0.48*y(k-4)+0.01*randn;
end

plot(y);
xlabel('Time');
ylabel('Amplitude');
grid on

% AR Markov model parameters
m=2;  % number of hidden states
p=2;  % AR model order

% Initial probabilities vector
c=rand(m,1);
c=c/sum(c);

% States transition matrix
A=rand(m,m);
A=A./(sum(A,1)'*ones(1,m));

% Initialize Sigma
sigma=1;

% Initialize model parameters ((p+1) x m) in a way as to 
% force model stability.
for i=1:m,
    poles=rand(1,p);
    eqx=poly(poles);
    alfa(:,i)=[rand -eqx(2:end)].';
end

% Run the Hamilton algorithm
[P1,P2]=hamilton_algorithm(A,c,y,alfa,sigma);

% Run the Kim algorithm
P3=kim_algorithm(A,P1,P2,p);