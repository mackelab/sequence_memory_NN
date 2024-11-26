function fig1d_stats(datapath)
% compare mean activity in differen task windows based on Zscored spiking
% activity
load(datapath,'bindat','convdat','sortinfo','expinfo');
% load times for averagaing windows 
ptimes=loadparams_avgtimes; % load times for averaging
tvec=((1:size(bindat,1)) - 1000)-1;% subtract baseline 1000ms to get time vector

npos=4; % number of stimulus positions
% mean sample activity
sample=[];
% mean delay
delay=[];
% mean probe
probe=[];

for stimpos=1:npos
    locdat=convdat(:,sortinfo(:,3)==stimpos);
    sample(:,stimpos)=squeeze(mean(locdat((ptimes.twins(stimpos,1):ptimes.twins(stimpos,2))+abs(tvec(1)),:)));
    delay(:,stimpos)=squeeze(mean(locdat((ptimes.twindelay(1):ptimes.twindelay(2))+abs(tvec(1)),:)));
    probe(:,stimpos)=squeeze(mean(locdat((ptimes.twinpanel(1):ptimes.twinpanel(2))+abs(tvec(1)),:)));
end

[pa_su(1),b,c]=anova1(sample,[],'off');
F_su(1)=b{2,5};
% print to screen
sprintf('Evoked Response F and p are %0.5g and %0.5g\n',F_su(1),pa_su(1))
