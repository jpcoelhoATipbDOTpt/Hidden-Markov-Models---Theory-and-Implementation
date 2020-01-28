function O=generate_data()

% Function used to generate training data from the HMM
% with parameters:
c=[0.6;0.4];
A=[0.7 0.3;0.4 0.6];
B=[0.1 0.4 0.5;0.7 0.2 0.1];

% Data generation considering route_nbr different routes 
% and observ_nbr observations per route
route_nbr=2;
observ_nbr=20;

% Initialize matrices:
H=zeros(route_nbr,observ_nbr); % Hidden
O=zeros(route_nbr,observ_nbr); % Observable

for route=1:route_nbr,
    H(route,1)=random_func(c);
    for p=2:observ_nbr,
        H(route,p)=random_func(A(H(route,p-1),:));
    end
    for p=1:observ_nbr,
        O(route,p)=random_func(B(H(route,p),:));
    end
end