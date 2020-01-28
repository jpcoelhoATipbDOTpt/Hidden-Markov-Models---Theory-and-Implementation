u1=[1;2];        % means
u2=[2;1];     
v1=[.1 0;0 .2];  % variances (positive definite)
v2=[.3 0;0 .2];
w=[0.4 0.6];     % weights
% Generate with this pdf a training vector:
N=500;           % Number of samples
M=2;
MaxIter=100;
x=generate_2D_MoG_Data(u1,v1,u2,v2,w,N);
%--------------------------------------------Training Algorithm
[sol_W,sol_Mu,sol_Sigma]=TrainMdG(x,M,MaxIter);
%--------------------------------------------------------------
v(1,:,:)=v1;v(2,:,:)=v2;
plotEM2(w,u,v,sol_W,sol_Mu,sol_Sigma);
