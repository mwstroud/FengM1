% close all; clear all
% xx=load('P2SOMEPSP');
% start=580;stop=850;
% dt=0.025;pp=[round(start/dt):round(stop/dt)];
% figure;plot(xx(pp,1),xx(pp,2),'black');
% plot(xx(pp,1),xx(pp,2),'black','LineWidth',2);
% box off

close all; clear all
xx=load('MC_CPcurrent_low');
% start=1;stop=850;
% dt=0.025;pp=[round(start/dt):round(stop/dt)];
figure;plot(xx(:,1),xx(:,2),'black');
plot(xx(pp,1),xx(pp,2),'black','LineWidth',2);
box off;