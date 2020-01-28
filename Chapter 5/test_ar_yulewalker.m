% Initialization
clear all;
clc;
close all;

N = 20000;                % Number of data samples
e = randn(1,N+200);       % Generate N->(0,1)
e = (e-mean(e)) / std(e); % Gaussian random var.
bootstrapsize = 200;      % Nbr. of bootstrap iter.

% Generate time series
lambda = zeros(1,N+200);  % Initialize lambda
alfa0=2.3;                % True parameters

% Poles location
alfa = -poly([0.95 0.5+0.5i 0.5-0.5i]); 

% Generate time series
for k=4:N+200,                       
    lambda(k) = alfa0+alfa(2)*lambda(k-1)+...
                alfa(3)*lambda(k-2)+...
                alfa(4)*lambda(k-3)+e(k);
end

% Remove initial transitory values
lambda(1:200)=[];                    
e(1:200)=[];              % (200 first samples)

% Parameters estimation by Yule-Walker method
alfae = aryulewalker(lambda,3);

% Parameters standard deviation estimation by bootstrap
estimation(1:3)=lambda(1:3); 
for k=4:N,
    % Observation estimation
    estimation(k)=alfae(1)+alfae(2)*lambda(k-1)+...
         alfae(3)*lambda(k-2)+alfae(4)*lambda(k-3);
    % One step ahead prediction error
    e(k)=lambda(k)-estimation(k);    
end
par(1,:)=alfae;           % Estimated parameters
for l=2:bootstrapsize
    inx=randperm(N);      % Randomize the error
    e=e(inx);
    % Compute new observations...     
    estimation(1:3)=lambda(1:3);     
    for k=4:length(lambda)
        estimation(k)=alfae(1) +...
         alfae(2)*estimation(k-1) +...
         alfae(3)*estimation(k-2) + ...
         alfae(4)*estimation(k-3) + e(k);
    end
    % Estimate new coefficients
    par(l,:) = aryulewalker(estimation,3); 
end

% Compute the parameters standard deviation
std(par)