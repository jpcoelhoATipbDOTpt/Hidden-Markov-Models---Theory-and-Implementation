% Used to test the Rosen algorithm
A=[1 1;1 5;-1 0;0 -1];
b=[2;5;0;0];
C=[];
d=[];
gradF='rosen_test_func';
x0=[-2;-2];
MaxIter=100;
x=rosen(gradF,x0,A,b,C,d,MaxIter);
% ---------------------------------------------------PLOTS
[x1,x2]=meshgrid(-1:.01:2);
y=double((x1+x2<=2)& (x1+5*x2<=5) & x1>=0 & x2>=0);
pcolor(x1,x2,y);  % Admissible region
shading interp;       
xlabel('x1');
ylabel('x2');
hold on;
z=2*x1.^2+2*x2.^2-2*x1.*x2-4*x1-6*x2; % Objective function
surfl(x1,x2,z);
shading interp;
% Solution
z=2*x(1).^2+2*x(2).^2-2*x(1)*x(2)-4*x(1)-6*x(2);
plot3(x(1),x(2),0,'kx')