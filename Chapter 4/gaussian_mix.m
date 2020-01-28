% Generate de Gaussian mixture
x=-5:.1:4;
y1=gaussian(x,0,1);         % Gaussian #1
y2=gaussian(x,1,0.5);       % Gaussian #2
y3=gaussian(x,2,0.1);       % Gaussian #3
w=[0.7 0.2 0.1];            % Relative weights

yt=w(1)*y1+w(2)*y2+w(3)*y3; % Gaussian mixture

% Generate 5000 samples from 3 distinct Gaussian
% according to a distribution function w
mix=random_func(w,5000);
h=zeros(1,5000);
for i=1:5000,
    if mix(i)==1,
        h(i)=0 + 1.*randn(1,1);
    elseif mix(i)==2,
        h(i)=1 + 0.5.*randn(1,1);
    else
        h(i)=2 + 0.1.*randn(1,1);
    end
end

% Present the data
plot(x,yt)
hold on
histogram(h,100,'Normalization','pdf')