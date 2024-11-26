function [params]=loadparams_ER()
% function [params]=loadparams_ER()


% p value 
params.crit_p=1e-3;
params.baseline=500;
% how long is window cut out for plotting
params.cutwindowlength=1000;
params.baseline_window=[500 1000];
params.binwidth=100;
params.s=100;
params.e=800;
params.windows(:,1)=[params.s:params.binwidth/2:params.e-100]';
params.windows(:,2)=[params.s+100:params.binwidth/2:params.e]';
params.mult=diff(params.baseline_window)/params.binwidth;
params.ncomp=size(params.windows,1);
params.nconsec=1;
params.numbinslat=3; % how many bins consecutively sign for latency estimation?
% params for latency analysis
params.binwidth_lat=18;
params.e_lat=700;% stop comparing 700ms after stim onset
params.windows_lat(:,1)=[params.s:params.binwidth_lat/2:params.e_lat-100]';
params.windows_lat(:,2)=[params.s+100:params.binwidth_lat/2:params.e_lat]';
% params for evaluating significance for PHC and rest
params.pwin=[170 550;250 600];

% average window for evoked response comparison among stimuli (kruskal
% wallis comparison between stimuli)
params.tavg=[200 500];


