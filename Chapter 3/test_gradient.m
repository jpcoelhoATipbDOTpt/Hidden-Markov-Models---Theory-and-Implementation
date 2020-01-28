function test_gradient(O)

% First considering only update in matrix A. B and c are known and 
% A is randomly initialized.
A=rand(2,2);
A=A./(sum(A,2)*ones(1,2));
B=[0.1 0.4 0.5;0.7 0.2 0.1];
c=[0.6;0.4];

% Learning rate
eta=0.005;

% Initialize likelihood vector. Each row is related to one of the
% 3 probability matrices/vectors
Likelihood=zeros(3,100);

% Iterations
for j=1:100
    [~,OMEGA]=gradientLA_norm(A,B,O,c); % Compute gradient
                                        % (normalized version)
    A=A+eta*OMEGA; % Update matrix A
    [~,lik]=forward_algorithm_norm(A,B,O,c); % Compute performance
    Likelihood(1,j)=lik;
end

% Second, considering only update in matrix B. A and c are known 
% and B is randomly initialized.
A=[0.7 0.3;0.4 0.6];
B=rand(2,3);
B=B./(sum(B,2)*ones(1,3));
c=[0.6;0.4];

% Iterations
for j=1:100
    [~,OMEGA]=gradientLB_norm(A,B,O,c); % Compute gradient 
                                        % (normalized version)
    B=B+eta*OMEGA; % Update matrix B
    [~,lik]=forward_algorithm_norm(A,B,O,c); % Compute performance
    Likelihood(2,j)=lik;
end

% Third, considering only update in vector c. A and B are known 
% and c is randomly initialized.
A=[0.7 0.3;0.4 0.6];
B=[0.1 0.4 0.5;0.7 0.2 0.1];
c=rand(2,1);
c=c./sum(c);

% Iterations
for j=1:100
    [dLdc]=gradientLc(B,O); % Compute gradient
    c=c+eta*dLdc; % Update vector c
    [~,lik]=forward_algorithm_norm(A,B,O,c); % Compute performance
    Likelihood(3,j)=lik;
end

% Final plots
plot(1:100,Likelihood(1,:));
xlabel('Iteration');
ylabel('Log-Likelihood');
grid on;
figure
plot(1:100,Likelihood(2,:));
xlabel('Iteration');
ylabel('Log-Likelihood');
grid on;
figure
plot(1:100,Likelihood(3,:));
xlabel('Iteration');
ylabel('Log-Likelihood');
grid on