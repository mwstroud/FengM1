function [ ] = DrawTheta(expSig, fftHan)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


areaCode = [0 1 2 3 4 5 6];
OutPath = '../../';
% OutPath = '..\..\Outputs2\explore\Gamma\';
tstop = 3000;
windowLen = 1024; %800/tstep;


% CA3
i = 1;
fileNames = dir([OutPath 'bmtk_' expSig '_SpikeTime' num2str(areaCode(i)) '.txt']);
fn = {fileNames.name}
cxxs = [];
pxxs = []; spikes = [];
% SpikeTimes = SpikeTimes( SpikeTimes(:, 1) == 1, :);
for sname = fn
    SpikeTimes = importdata([OutPath sname{1}]);

    ev = zeros(1, tstop);
    rSPT = round(SpikeTimes(:,2));
    for ST = 1:length(rSPT)
       ev(rSPT(ST)) = ev(rSPT(ST)) + 1; 
    end
 
        nID = unique(SpikeTimes(:,1));
    spikes(end+1) = length(SpikeTimes(:,2));
    [f,Pxxn,tvect,Cxx] = psautospk(ev, 1, windowLen, bartlett(windowLen), windowLen/2) ;
 
    pxxs(:,end+1) = Pxxn;
    cxxs(:, end+1)= Cxx;
 
end
 
  [mean(spikes) std(spikes)]/(63*30)

xx = 0:size(pxxs,1)-1;
yy = mean(pxxs,2);
sy = std(pxxs');
 
% hsubplots(FFTavg) = subplot (yplots,xplots, FFTavg); % second subplot
axes(fftHan)
rr = 1000;
boundedline(xx', yy/rr, sy/rr, 'k')

% chil = get(hsubplots(fftHan), 'children');
% set (chil(2), ...
%     'Color'         , co , ...
%     'LineWidth'     , 1             );

        xlim([0  50]);
         set(gca, 'XTick'       , 0:10:50);
                  set(gca, 'XTick'       , 0:10:50);  
         set(gca, 'XTickLabel',['  '; '10' ; '  '; '30'; '  '; '50'])
             



 
end