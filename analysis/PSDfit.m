clear all;close all;clc;
fs=1000;nfft=1024/2;
path='';strnum='33';str1='LFP_elec_0_',str2='_shortburst';str3='';

LFP=load(fullfile(path,strcat(str1,strnum,str2,str3)));
tstop=length(LFP);
burst_len=200;rest_len=500;
burst_start=[1:burst_len+rest_len:tstop];
burst_stop=[1+burst_len:burst_len+rest_len:tstop];

burstflag=3; %%%%1 means burst period. 2 means resting peroid, 3 means all

lfp_burst=[];lfp_rest=[];
if burstflag==1;
    for i=1:numel(burst_stop)
    lfp_burst=[lfp_burst;LFP(burst_start(i):burst_stop(i))];
    end
    lfp_test=lfp_burst;
elseif burstflag==2;
     for i=1:numel(burst_stop)-1
    lfp_rest=[lfp_rest;LFP(burst_stop(i):burst_start(i+1))];
     end
     lfp_test=lfp_rest;
elseif burstflag==3;
     lfp_test=LFP;
end

LFP=lfp_test;
LFP=LFP(1:end,1);
[b,a] = butter(2,1.17/500,'high');
% 
% LFP=filtfilt(b,a,LFP);

% matlabWork='';
% freqs = 2.^(0:0.25:8);
% powerFunc = @(x)mean(abs(x.Spectrum),2);
% specPoiss = cellfun(@(x)powerFunc(WaveletSpec((filtfilt(b,a,x)*1.8*100)/1000, ...
%   0.001, freqs, 'NORMALIZE')),mat2cell(LFP,length(LFP),1),'uniformoutput',false);
% pxx=(specPoiss{1,1});f=freqs;
[pxx,f] = pwelch(LFP,nfft,0,nfft,fs);
figure
plot(f,pxx)
 set(gca, 'YScale', 'log');
 set(gca, 'XScale', 'log');
fit_range = [1,90];  % [10,200]
f_fit = find(f>=fit_range(1)&f<=fit_range(2));
% use evenly distributed samples in log scale to fit
f_log = log(f(f_fit));
% f_log = (f(f_fit));

f_pts = exp(linspace(f_log(1),f_log(end),200));
[~,f_fit] = min(abs(bsxfun(@minus,f_pts,f)),[],1);
f_fit = unique(f_fit);

outlier_threshold = 1;    % 1. outlier MAD threshold, multiple of 1 sigma.
[fit_a,fit_b,outliers] = powerfit(f(f_fit),pxx(f_fit),1,outlier_threshold,1);
outliers = f_fit(outliers);
pxxf = fit_a*f.^fit_b;

figure; hold on;
plot(f,pow2db(pxx),'b');
plot(f,pow2db(pxxf),'r--','LineWidth',2);
plot(f(outliers),pxx(outliers),'m.');
set(gca,'xscale','log');

figure(311)
title('1/f residue')
plot(f,pow2db(pxx)-pow2db(pxxf),'blue','LineWidth',2);
set(gca,'xscale','log');xlim([5,500]);grid on;ylim([-10,7]);
y = 0;
%line([1,500],[y,y],'Color','red','LineStyle','--')


%%%%%%%%%%%%%%%%%%%%%%%%%

clear LFP;
str2='_longburst';
LFP=load(fullfile(path,strcat(str1,strnum,str2,str3)));

tstop=length(LFP);
burst_len=1000;rest_len=500;
burst_start=[1:burst_len+rest_len:tstop];
burst_stop=[1+burst_len:burst_len+rest_len:tstop];


lfp_burst=[];lfp_rest=[];
if burstflag==1;
    for i=1:numel(burst_stop)
    lfp_burst=[lfp_burst;LFP(burst_start(i):burst_stop(i))];
    end
    lfp_test=lfp_burst;
elseif burstflag==2;
     for i=1:numel(burst_stop)-1
    lfp_rest=[lfp_rest;LFP(burst_stop(i):burst_start(i+1))];
     end
     lfp_test=lfp_rest;
elseif burstflag==3;
     lfp_test=LFP;
end

LFP=lfp_test;
LFP=LFP(1:end,1);
[b,a] = butter(2,1.17/500,'high');
% 
% LFP=filtfilt(b,a,LFP);

% matlabWork='';
% freqs = 2.^(0:0.25:8);
% powerFunc = @(x)mean(abs(x.Spectrum),2);
% specPoiss = cellfun(@(x)powerFunc(WaveletSpec((filtfilt(b,a,x)*1.8*100)/1000, ...
%   0.001, freqs, 'NORMALIZE')),mat2cell(LFP,length(LFP),1),'uniformoutput',false);
% pxx=(specPoiss{1,1});f=freqs;
[pxx,f] = pwelch(LFP,nfft,0,nfft,fs);
figure
plot(f,pxx)
 set(gca, 'YScale', 'log');
 set(gca, 'XScale', 'log');
fit_range = [1,90];  % [10,200]
f_fit = find(f>=fit_range(1)&f<=fit_range(2));
% use evenly distributed samples in log scale to fit
f_log = log(f(f_fit));
% f_log = (f(f_fit));

f_pts = exp(linspace(f_log(1),f_log(end),200));
[~,f_fit] = min(abs(bsxfun(@minus,f_pts,f)),[],1);
f_fit = unique(f_fit);

outlier_threshold = 1;    % 1. outlier MAD threshold, multiple of 1 sigma.
[fit_a,fit_b,outliers] = powerfit(f(f_fit),pxx(f_fit),1,outlier_threshold,1);
outliers = f_fit(outliers);
pxxf = fit_a*f.^fit_b;

figure; hold on;
plot(f,pow2db(pxx),'b');
plot(f,pow2db(pxxf),'r--','LineWidth',2);
plot(f(outliers),pxx(outliers),'m.');
set(gca,'xscale','log');

figure(311)
title('1/f residue')
hold on;plot(f,pow2db(pxx)-pow2db(pxxf),'red','LineWidth',2);
set(gca,'xscale','log');xlim([5,500]);grid on;ylim([-10,4]);legend('short burst input','long burst input')
y = 0;
line([1,500],[y,y],'Color','black','LineStyle','--','LineWidth',2)


