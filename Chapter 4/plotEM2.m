function plotEM2(W,Mu,Sigma,sol_W,sol_Mu,sol_Sigma)
[~,pdf]=drawGauss2(W,Mu,Sigma);
figure;
[xx,pdf2]=drawGauss2(sol_W(1,:),sol_Mu(1,:,:),sol_Sigma(1,:,:,:));
N1=sqrt(length(xx));
x1=reshape(xx(1,:),N1,N1);
x2=reshape(xx(2,:),N1,N1);
pdf=reshape(pdf,N1,N1);
pdf2=reshape(pdf2,N1,N1);
h=surfl(x1,x2,pdf);
set(h,'FaceAlpha',0.4);
shading interp
hold on;
h=surfl(x1,x2,pdf2);
shading interp
colormap gray;
set(h,'EdgeColor',[0.5 0.5 0.5]);
axis tight
%-----------------------------------------------------Sub Function
[xx,pdf2]=
drawGauss2(sol_W(end,:),sol_Mu(end,:,:),sol_Sigma(end,:,:,:));
figure;
N1=sqrt(length(xx));
x1=reshape(xx(1,:),N1,N1);
x2=reshape(xx(2,:),N1,N1);
pdf=reshape(pdf,N1,N1);
pdf2=reshape(pdf2,N1,N1);
h=surfl(x1,x2,pdf);
set(h,'FaceAlpha',0.4);
shading interp
hold on;
h=surfl(x1,x2,pdf2);
shading interp
colormap gray;
set(h,'EdgeColor',[0.5 0.5 0.5]);
axis tight
