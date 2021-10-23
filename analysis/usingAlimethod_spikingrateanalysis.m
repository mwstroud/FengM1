clear all;close all;
clc;
addpath('F:\Feng_gammNSG')
data = load('data_32_longburst');
% data=data(data(:,1)<=5000,:);
tstop=max(data(:,1));
windowLen = 1024;
ev_I = zeros(1, tstop);
        
figure (121212)
upscale=1;
num_P=800*upscale;inter_num=100*upscale;som_num=100*upscale;
all_num=num_P+inter_num+som_num;
        
        
data_P=data(find(data(:,2)<num_P),:);
data_I=data(find(data(:,2)>=num_P&data(:,2)<num_P+inter_num),:);
data_SOM=data(find(data(:,2)>=num_P+inter_num&data(:,2)<all_num),:);

plot(data_P(:,1),data_P(:,2)+1,'blue.');
hold on; plot(data_I(:,1),data_I(:,2)+1,'r.');  
hold on; plot(data_SOM(:,1),data_SOM(:,2)+1,'g.'); 

data_P_analysis=data_P;
data_I_analysis=data_I;
data_SOM_analysis=data_SOM;

exclusive_time=1; % (ms)Excluding time period for analysis


rSPT = round(data_I_analysis(:,1)*1);
rSPT_sort=sort(rSPT);
rSPT_unique=unique(rSPT_sort);
[n, bin] = histc(rSPT_sort, rSPT_unique);
ev_I=zeros(tstop,1);
ev_I(rSPT_unique)=n;

%     for ST = 1:length(rSPT)
%        ev_I(rSPT(ST)) = ev_I(rSPT(ST)) + 1; 
%     end
   
%  Excluding the first 300 ms from calculation
 ev_I(1:exclusive_time) = [];
     
 [f,Pxxn,tvect,Cxx] = psautospk(ev_I, 1, windowLen, bartlett(windowLen), windowLen/2, 'none') ;
 figure (1212)
 subplot(3,1,2)  
 plot(f,Pxxn,'r', 'LineWidth', 2);
 set(gca,'yscale','linear'); legend('PSD for ITN firing rate');xlabel('Hz');ylabel('PSD');
 
    
rSPT = round(data_P_analysis(:,1)*1);

rSPT_sort=sort(rSPT);
rSPT_unique=unique(rSPT_sort);
[n, bin] = histc(rSPT_sort, rSPT_unique);
ev_P=zeros(tstop,1);
ev_P(rSPT_unique)=n;

     ev_P(1:exclusive_time) = [];
     
    [f,Pxxn,tvect,Cxx] = psautospk(ev_P, 1, windowLen, bartlett(windowLen), windowLen/2, 'none') ;

     subplot(3,1,3)  
     plot(f,Pxxn,'blue', 'LineWidth', 2);
     set(gca,'yscale','linear'); legend('PSD for PN firing rate');xlabel('Hz');ylabel('PSD');
     
    

rSPT = round(data_SOM_analysis(:,1)*1);

rSPT_sort=sort(rSPT);
rSPT_unique=unique(rSPT_sort);
[n, bin] = histc(rSPT_sort, rSPT_unique);
ev=zeros(tstop,1);
ev(rSPT_unique)=n;
     ev(1:exclusive_time) = [];
     
    [f,Pxxn,tvect,Cxx] = psautospk(ev, 1, windowLen, bartlett(windowLen), windowLen/2, 'none') ;

     subplot(3,1,1)  
     plot(f,Pxxn,'green', 'LineWidth', 2);
     set(gca,'yscale','linear'); legend('PSD for SOM firing rate');xlabel('Hz');ylabel('PSD');
     
     
%%%%plot moving-windowed firing rate     
     
binSize = 20; %%ms
spiketimes = data_I_analysis(:, 1); 
timeMax = max(spiketimes);
timeMax = timeMax - binSize;
binSpikeCount = [];
tt_up=[0:1:tstop]+(binSize/2);

N_up = histc(spiketimes, tt_up);

movingSum = conv(N_up, ones(1, binSize));% N_up = histc(spiketimes, tt_up);

figure (1212123)
hold on; plot(movingSum*1000/binSize,'r');
axis tight

spiketimes = data_P_analysis(:, 1); 
% timeMax = max(spiketimes);
% timeMax = timeMax - binSize;
% binSpikeCount = [];
N_up = histc(spiketimes, tt_up);

movingSum = conv(N_up, ones(1, binSize));% N_up = histc(spiketimes, tt_up);

hold on; plot(movingSum*1000/binSize,'b');% plot(firingRate) 


spiketimes = data_SOM_analysis(:, 1); 
% timeMax = max(spiketimes);
% timeMax = timeMax - binSize;
% binSpikeCount = [];
N_up = histc(spiketimes, tt_up);

movingSum = conv(N_up, ones(1, binSize));% N_up = histc(spiketimes, tt_up);

hold on; plot(movingSum*1000/binSize,'g');% plot(firingRate) 

legend('ITN firing rate by moving window','PNs firing rate by moving window','SOMs firing rate by moving window')
axis tight;


%%%%plot binned firing rate     
% figure (300)
% %  Excluding the first 300 ms from calculation
% %  firingRate(1:exclusive_time)=[];
% 
% subplot(3,1,1)
% firingRate_nonoverlap=ev_I;
% bar(firingRate_nonoverlap,'red'); axis tight; 
% %set(gca,'XTick',[exclusive_time:5:tstop]);
% hold on;
% firingRate_nonoverlap=ev_P;
% bar(firingRate_nonoverlap,'blue'); axis tight;
% legend('ITN spiking rate in 1ms bin','PN spiking rate in 1ms bin');
% 
% 
% subplot(3,1,2)
% firingRate_nonoverlap=ev_I;
% bar(firingRate_nonoverlap,'red'); axis tight;
% legend('ITN spiking rate in 1ms bin')
% 
% subplot(3,1,3)
% firingRate_nonoverlap=ev_P;
% bar(firingRate_nonoverlap,'blue'); axis tight;
% legend('PN spiking rate in 1ms bin')
% 




% title('Spiking rate(by sliding window)/(Hz)');



