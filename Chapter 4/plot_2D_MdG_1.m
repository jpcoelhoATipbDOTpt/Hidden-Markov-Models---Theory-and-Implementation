dx = .1;                       % Space resolution
x  = -1:dx:4;                  % Space limits
N = [length(x) length(x)];     % Number of points
[x1,x2]=meshgrid(x);           % Create mesh
xx=[x1(:) x2(:)]';             % 2D data for plot
Nt = length(xx);               % Number of points
u1=[1;2];                      % Mean for Gauss. 1
v1=[.1 0;0 .2];                % Cov. for Gauss. 1
u2=[2;1];                      % Mean for Gauss. 2
v2=[.3 0;0 .2];                % Cov. for Gauss. 2
w =[0.4 0.6];                  % Mixing weights
pdf = zeros(1,Nt);             % Mixed pdf
%-------------------------------------------------
for j=1:Nt;
  pdf(j)=w(1)/(sqrt(det(2*pi*v1)))*...
  exp(-0.5*(((xx(:,j)-u1)'/v1)*(xx(:,j)-u1)))+...
  w(2)/(sqrt(det(2*pi*v2)))*...
  exp(-0.5*(((xx(:,j)-u2)'/v2)*(xx(:,j)-u2)));
end
%-------------------------------------------------
x1=reshape(xx(1,:),N);         
x2=reshape(xx(2,:),N);         
pdf=reshape(pdf,N);            
%-------------------------------------------------
surfl(x1,x2,pdf)               % Plot PDF
shading interp                 %
colormap gray                  %
axis tight                     %
xlabel('x_1');ylabel('x_2');   %