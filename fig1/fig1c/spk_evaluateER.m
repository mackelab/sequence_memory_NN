function [ERstats,params,outdat]=spk_evaluateER(data,expinfo,params,clust_info)
% this plot evaluates the evoked response after stimulus onset across all
% stimuli ids shown. It uses predefined parameters from loadparams_ER.m
% that outputs params variable
% function [ERstats,params,outdat]=spk_evaluateER(data,expinfo,params_ER,clust_info)

% Experimental Parameters
nstims=numel(expinfo.stimuli.stimnum);
stimids=expinfo.stimuli.stimnum;
stimpos=1:size(expinfo.times.stimon,2);
npos=numel(stimpos);

% cut out evoked response window
% used later for plotting
% 4 diffeclrent times for each position
cutstart=median(expinfo.times.stimon);
cutstop=median(expinfo.times.stimon)+params.cutwindowlength;

% timewindow to evaluate p-values is different for regions:

switch clust_info.sitename{:}
    case 'LPHC'
        pwin=params.pwin(1,:);
    otherwise
        pwin=params.pwin(2,:);
end


% matrix with info about trials
trialinfo=[];

P_sr=zeros(nstims,size(params.windows,1));
Tstamps=cell(nstims,1);
outdat=struct('trialsum',[],'baseline',[],'chron_dat',[]);
ERstats=struct('P_sr',[]);

for i=1:nstims
    %
    basedat=[];
    stimdat=[];
    sortinfo=[];
    
    for j=1:npos
        % sortdata according to stim id, j'th position
        [sortdata,tempsort]=sorttrials(data,expinfo, stimids(i),stimpos(j));
        
        % cut baseline data
        [base]=spk_cutdata(sortdata,params.baseline_window);
        
        % cut evoked response 0 to 1000ms post stim_onset at each position
        [cdata]=spk_cutdata(sortdata,[cutstart(j) cutstop(j)]);
        
        % bin data per trial
        [bindat]=spk_timestobin(cdata,[cutstart(j)+1 cutstop(j)]);
        [bindatbase]=spk_timestobin(base,[params.baseline_window(1)+1 params.baseline_window(2)]);
        
        
        % collect binned data
        stimdat=[stimdat;bindat];
        
        basedat=[basedat;bindatbase];
        
        sortinfo=[sortinfo;tempsort];
        
        % save timestamps for histogram plotting later
        % concatenate baselinedata und stimulus evoked data timestamps in
        % one big vector
        % align timestamps of evoked response window so that 0 corresponds to stimulus onset
        % from baseline only
        Tstamps{i}=[Tstamps{i};[cell2mat(base)-cutstart(1) ;...
            cell2mat(cellfun(@(cdata) cdata-cutstart(j),cdata,'UniformOutput',false))]];
        
        
    end % end positions
    
    
    % put trials again to original chronological order
    Ntrials=size(sortinfo,1);
    [s,ix]=sort(sortinfo(:,1));
    chron_dat=[zeros(size(basedat)) zeros(size(stimdat))];
    
    % time stamps 1ms bins all trials per stimuli
    chron_dat=[basedat(ix,:) stimdat(ix,:)];
    
    % sort trialinfo according to sortindex
    % first column is trial index, second column is stimulus id, third is,
    % forth is correct or not
    
    % stimulus position
    trialinfo=[trialinfo;[sortinfo(ix,:) expinfo.resp.correct(ix)]];
    % if you want reverse oder in plotting need to swap matrix here...
    
    % sum of spikes in baseline window, convert to Hz
    basesum=sum(basedat,2)*(1000/diff(params.baseline_window));
    p_sr=zeros(1,length(params.windows));
    
    for k=1:length(params.windows)
        
        % count spikes in response windows, convert to Hz (spk/s)
        trialsum(:,k)=sum(stimdat(:,params.windows(k,1):params.windows(k,2)),2).*(1000/params.binwidth);
        
        % for each window
        % do signed rank test on baseline vs. stim evoked data
        [p_sr(k),temph]=signrank(trialsum(:,k),basesum,'tail','right');
        
    end
    
    for k=1:length(params.windows_lat)
        
        % count spikes in response windows, convert to Hz (spk/s)
        trialsum_lat(:,k)=sum(stimdat(:,params.windows_lat(k,1):params.windows_lat(k,2)),2).*(1000/params.binwidth_lat);
        
        % for each window
        % do signed rank test on baseline vs. stim evoked data
        [p_sr_lat(i,k),h_lat(i,k)]=signrank(trialsum_lat(:,k),basesum,'tail','right');
        
        
    end
    
    % correct p values Simes
    
    
    % latency estimation
    % look for large bin with least p value / sign response
    
    % look for first three consecutive bins with sign p value
    
    
    % correct p values by number of comparisons
    % sort p-vals
    [sp,idx]=sort(p_sr);
    % correct them w/ N of comps and put them back at bin
    % they belong
    p_sr(idx)=(params.ncomp./(1:params.ncomp)).*sp; % Simes correction for multiple comparison
    P_sr(i,:)=p_sr;
    
    
    
    % save binned data per stimulus (1ms bins)
    
    % save trial-summed data per bin
    outdat.trialsum(:,:,i)=trialsum;
    outdat.baseline(:,:,i)=basesum;
    % both baseline and stimulus response (artifically concatenated next to
    % eachother)
    outdat.chron_dat(:,:,i)=chron_dat;
    outdat.trialinfo(:,:,i)=trialinfo;
    trialinfo=[];
    
end % end stimuli

ERstats.P_sr=P_sr;
outdat.Tstamps=Tstamps;


% find minimal p value within given time windows
% index into data
% tempix=find(params.windows(:,1)==pwin(1)):find(params.windows(:,2)==pwin(2));
[m,tempix]=min(abs([params.windows(:,1)-pwin(1) params.windows(:,2)-pwin(2)]));
tempix=tempix(1):tempix(2);


% tempwindows to evaluate
tempwin=params.windows(find(params.windows(:,1)==pwin(1)):find(params.windows(:,2)==pwin(2)),:);
tempwin=params.windows(tempix,:);

% minimum p value within time frame
[pER,ix]=min(ERstats.P_sr(:,tempix),[],2);
% window with minimum p-value (only important for plotting)
pwin=tempwin(ix,:);

% 
% tempix=find(params.windows(:,1)==pwin(1)):find(params.windows(:,2)==pwin(2));
% % tempwindows to evaluate
% tempwin=params.windows(find(params.windows(:,1)==pwin(1)):find(params.windows(:,2)==pwin(2)),:);
% % minimum p value within time frame
% [pER,ix]=min(ERstats.P_sr(:,tempix),[],2);
% % window with minimum p-value (only important for plotting)
% pwin=tempwin(ix,:);

% stimulus id with highest response increase 
beststim=find(pER==min(pER));




% add to ERstats
ERstats.pER=pER;
ERstats.pwin=pwin;
ERstats.beststim=beststim;
ERstats.p_sr_lat=p_sr_lat;
ERstats.h_lat=h_lat;
ERstats.beststim=beststim;


% assess latency of response to best stimulus (lowest p)
%window_latency=ERstats.pwin(ERstats.pER==min(ERstats.pER),:);
%[x,ixx]=min(abs(params.windows_lat(:,1)-window_latency(1)));
% find bin for which first sign resp

%%






