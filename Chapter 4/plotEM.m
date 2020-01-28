function plotEM(W,Mu,Sigma,sol_W,sol_Mu,sol_Sigma)
[~,pdf1]=drawGauss(W,Mu,Sigma);
figure
subplot(2,1,1);
[xx,pdf2]=drawGauss(sol_W(1,:),sol_Mu(1,:),sol_Sigma(1,:));
plot(xx,pdf1,xx,pdf2); axis([-1 5 0 0.8]);xlabel('x');
legend('real distribution','estimated distribution'); 
text(-0.8,0.7,'Iteration #1');
subplot(2,1,2);
[xx,pdf2]=drawGauss(sol_W(end,:),sol_Mu(end,:),sol_Sigma(end,:));
plot(xx,pdf1,xx,pdf2);axis([-1 5 0 0.8]);xlabel('x');
legend('real distribution','estimated distribution');
text(-0.8,0.7,'Iteration #1000');
