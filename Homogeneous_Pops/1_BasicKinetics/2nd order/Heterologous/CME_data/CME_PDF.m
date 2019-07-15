% Calculate discretely-sampled PDF of transition events based on mean CME trajectory.
% clear;

Ao = 50;
Bo = ceil(0.7*Ao);               % B is limiting!  

load(['CME_A=' num2str(Ao) '_B=' num2str(Bo) '.mat']);

dt = 1/100; 
t_max = 10;    
bins = 20;

time_PDF = 0:t_max/bins:t_max;
%%
avg_B_PDF_temp = [Bo E];

time_temp = time_PDF / dt + 1;
time_temp(end) = t_max/dt;

avg_B_PDF = avg_B_PDF_temp(time_temp); 

dB_CME = abs(diff(avg_B_PDF)) / Bo;                     

time_PDFm = t_max/bins/2:t_max/bins:t_max;

% plot(time_PDFm,dB_CME','Color',[0 0.75 0],'LineWidth',2);

%% Calculate R^2 
CoefDet(Num1,dB_CME)