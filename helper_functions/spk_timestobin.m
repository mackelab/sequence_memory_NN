function [bindat]=spk_timestobin(trialdat,twindow,binsize)
% [bindat]=spk_timestobin(trialdat,twindow,binsize)
% converts spktimes of cellarray (spktimes per trials) into one matrix with 
% equal timebins, useful for rasterplotting
% input: timestamps of one unit, twindow is first and last entry of edges
% between which elements are counted
% if data aligned to fixation onset (default): edges are ms post fix_onset
% eg [100 2000] takes all spikes 100ms post fix_on until 2000ms post fix on
% if data realigned e.g. to stim_onset, 0 in tstamps corresponds to stim_onset
% need to adjust edges i.e. twindow - timeshift(ms)
% output: datamatrix with (ntrials, timebins) in 1ms resolution
% default: binsize=1ms;
% SL 

% dummy test
% v=zeros(100,1);
% v(1:20)=1;
% fv=find(v);
% hv=spk_timestobin({fv},[1 100]);
if nargin<3
binsize=1;
end
ntrials = size(trialdat,1);
bindat=zeros(ntrials,diff(twindow)+1);


for i=1:ntrials
    tempdat=histc(trialdat{i},twindow(1):binsize:twindow(2));
    if isempty(tempdat)
        tempdat=zeros(1,diff(twindow)+1);
    end
    bindat(i,:)=tempdat;
end

