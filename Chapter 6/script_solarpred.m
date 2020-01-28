clear all;
close all;
clc;

% Load matrix with 1440 x 361 of measured radiation
load('solardata.mat'); 

% Parameters Estimation data
Lambda = [sol.day(:,150);sol.day(:,151);sol.day(:,152)].'; 

% Constants
MaxIter=100; % maximum number of training iterations
m=3;         % number of hidden states
p=4;         % AR model order
h=10;        % prediction horizon

% Training
[A,c,alfa,sigma2]=ARMM_training(Lambda,m,p,MaxIter);

% Prediction
Lambda = solar.day(:,153).'; % Validation Data
[Lambda_pred,Lambda_meas] = ...
   ARMM_kstepahead(Lambda,A,c,alfa,sigma2,h);
   
% Results
figure(1);
plot([Lambda_meas(:,1) Lambda_pred(:,1)])
axis([0 1440 0 350])
xlabel('Time/min');
ylabel('Solar Irradiation/Wm^{-2}');
legend('Measured','one step ahead prediction');
error=Lambda_meas(:,1)-Lambda_pred(:,1);

% Compute prediction error
seq1=sqrt((error'*error)/length(error));  
figure(2);
plot([Lambda_meas(:,h) Lambda_pred(:,h)])
axis([0 1440 0 400])
xlabel('Time/min');
ylabel('Solar Irradiation/Wm^{-2}');
legend('Measured',[num2str(h) ' step ahead prediction']);
y=Lambda_meas(:,h);ye=Lambda_pred(:,h);
error=Lambda_meas(:,h)-Lambda_pred(:,h);

% Compute prediction error
seq60=sqrt((error'*error)/length(error));