function pdf=Gfit2D(x,W,Mu,Sigma)

L=size(Sigma);
S1=reshape(Sigma(1,:,:),L(2),L(3));
S2=reshape(Sigma(2,:,:),L(2),L(3));
M1=Mu(:,1);
M2=Mu(:,2);
for j=1:length(x);
    pdf=W(1)/(sqrt(det(2*pi*S1)))*...
    exp(-0.5*(((x(:,j)-M1)'/S1)*(x(:,j)-M1)))+...
    W(2)/(sqrt(det(2*pi*S2)))*...
    exp(-0.5*(((x(:,j)-M2)'/S2)*(x(:,j)-M2)));
end