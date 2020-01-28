function [sol_W,sol_Mu,sol_Sigma,Fit]=test_trainMdG2D
% Problem 2D @ 2MoG

% (c)2010 João Paulo Coelho
%PDF
dx=.1;
[x1,x2]=meshgrid(-1:dx:4);
xx=[x1(:) x2(:)]';
% Lines - define dimensions
% Columns - define Gaussian
u=[1 2;2 1]; % means
v1=[.1 0;0 .2];   % variances (positive definite)
v2=[.3 0;0 .2];
% weights
d=[0.4 0.6];
% -----------------PDF
for k=1:length(xx);
    pdf(k)=d(1)/(sqrt(det(2*pi*v1)))*exp(-0.5*(((xx(:,k)-u(:,1))'/v1)*(xx(:,k)-u(:,1))))+...
        d(2)/(sqrt(det(2*pi*v2)))*exp(-0.5*(((xx(:,k)-u(:,2))'/v2)*(xx(:,k)-u(:,2))));
end
N1=size(x1);
x1=reshape(xx(1,:),N1);
x2=reshape(xx(2,:),N1);
pdf=reshape(pdf,N1);

figure
surfl(x1,x2,pdf)
shading interp
colormap gray
axis tight
% view(-9,18)
xlabel('x_1');ylabel('x_2');

% Generate with this pdf a training vector:
N=500;
x=zeros(2,N);
for k=1:N,
    inx=rand;
    if inx<d(1);
        x(:,k)=u(:,1)+sqrt(v1)*randn(2,1);
    else
        x(:,k)=u(:,2)+sqrt(v2)*randn(2,1);
    end
end
%------------------------------------------Shows the data distribution
hold on
plot3(x(1,:),x(2,:),zeros(1,N),'.')
%-------------------------------------------------------Training Algorithm
[sol_W,sol_Mu,sol_Sigma]=TrainMdG(x,2,100);
%--------------------------------------------------------------------------

v(1,:,:)=v1;
v(2,:,:)=v2;
plotEM(d,u,v,sol_W,sol_Mu,sol_Sigma)
%----------------------------------------------------------------Sub Function
function plotEM(W,Mu,Sigma,sol_W,sol_Mu,sol_Sigma)
% (c)2010 João Paulo Coelho
[xx,pdf]=drawGauss(W,Mu,Sigma);
figure
% subplot(2,1,1);
[xx,pdf2]=drawGauss(sol_W(1,:),sol_Mu(1,:,:),sol_Sigma(1,:,:,:));
N1=sqrt(length(xx));
x1=reshape(xx(1,:),N1,N1);
x2=reshape(xx(2,:),N1,N1);
pdf=reshape(pdf,N1,N1);
pdf2=reshape(pdf2,N1,N1);
h=surfl(x1,x2,pdf);
set(h,'FaceAlpha',0.4)
shading interp
hold on
h=surfl(x1,x2,pdf2);
shading interp
colormap gray
set(h,'EdgeColor',[0.5 0.5 0.5])
axis tight

%----------------------------------------------------------------Sub Function
[xx,pdf2]=drawGauss(sol_W(end,:),sol_Mu(end,:,:),sol_Sigma(end,:,:,:));
% (c)2010 João Paulo Coelho
figure
% subplot(2,1,2);
N1=sqrt(length(xx));
x1=reshape(xx(1,:),N1,N1);
x2=reshape(xx(2,:),N1,N1);
pdf=reshape(pdf,N1,N1);
pdf2=reshape(pdf2,N1,N1);
h=surfl(x1,x2,pdf);
set(h,'FaceAlpha',0.4)
shading interp
hold on
h=surfl(x1,x2,pdf2);
shading interp
colormap gray
set(h,'EdgeColor',[0.5 0.5 0.5])
axis tight

%----------------------------------------------------------------Sub Function
function [xx,pdf]=drawGauss(W,Mu,Sigma)
% (c)2010 João Paulo Coelho
L=size(Sigma);
if length(L)>3
    v1=reshape(Sigma(1,1,:,:),L(3),L(4));
    v2=reshape(Sigma(1,2,:,:),L(3),L(4));
else
    v1=reshape(Sigma(1,:,:),L(2),L(3));
    v2=reshape(Sigma(2,:,:),L(2),L(3));
end
L=size(Mu);
if length(L)>2
    Mu=reshape(Mu(1,:,:),L(2),L(3));
end
dx=.1;
[x1,x2]=meshgrid(-1:dx:4);
xx=[x1(:) x2(:)]';
for k=1:length(xx);
    pdf(k)=W(1)/(sqrt(det(2*pi*v1)))*exp(-0.5*(((xx(:,k)-Mu(:,1))'/v1)*(xx(:,k)-Mu(:,1))))+...
        W(2)/(sqrt(det(2*pi*v2)))*exp(-0.5*(((xx(:,k)-Mu(:,2))'/v2)*(xx(:,k)-Mu(:,2))));
end
