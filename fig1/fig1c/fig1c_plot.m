function fig1c_plot(datafilepath,imagepath)
% provide path to data and images
load(datafilepath,'spkdat','clust_info','outdat','expinfo','params','ERstats', 'spikes', 'idx','unit');
nstimuli=numel(expinfo.stimuli.stimnum);

% parameters for plotting
plotparams.binsize=50;
plotparams.edgv=[-params.baseline:plotparams.binsize:params.cutwindowlength];
params.plotwaveform=1;
params.normHz{1}=1;params.normHz{2}=112;% ntrials per stimulus shown 

% image
% raster
% PSTH
% define positions to plot
if params.plotwaveform
    rasterpos{1}=repmat([0.03 0.23 0.43 0.63],1,2);
    rasterpos{2}=repmat([0.62 0.12],4,1);
    rasterpos{3}=repmat([0.17 0.28],8,1);
    
    histpos{1}=repmat([0.03 0.23 0.43 0.63],1,2);
    histpos{2}=repmat([0.53 0.03],4,1);
    histpos{3}=repmat([0.17 0.08],8,1);
    
    stimpos{1}=repmat([0.08 0.28 0.48 0.68],1,2);
    stimpos{2}=repmat([0.92 0.42],4,1);
    stimpos{3}=repmat([0.07 0.07],8,1);
    
    
else
    rasterpos{1}=repmat([0.03 0.27 0.52 0.77],1,2);
    rasterpos{2}=repmat([0.62 0.12],4,1);
    rasterpos{3}=repmat([0.17 0.28],8,1);
    
    histpos{1}=repmat([0.03 0.27 0.52 0.77],1,2);
    histpos{2}=repmat([0.53 0.03],4,1);
    histpos{3}=repmat([0.17 0.08],8,1);
    
    stimpos{1}=repmat([0.08 0.34 0.58 0.83],1,2);
    stimpos{2}=repmat([0.92 0.42],4,1);
    stimpos{3}=repmat([0.07 0.07],8,1);
end
histhandles=cell(nstimuli,1);
%
fh=figure;
set(gcf,'position',[590 35 1070 670]);
% loop across images
for s=1:nstimuli
    
    % stimulus image
    % load image
    tempimage=imread([imagepath expinfo.stimuli.stimnames{s}]);
    % plot stimulus image
    ah_s(s)=axes('position',[stimpos{1}(s) stimpos{2}(s) stimpos{3}(s,:)]);
    imshow(tempimage);
    
    
    % Raster plots
    chron_dat=squeeze(outdat.chron_dat(:,:,s));
    ah_r(s)=axes('position',[rasterpos{1}(s) rasterpos{2}(s) rasterpos{3}(s,:)]);
    plotRasters(chron_dat,'x',[1 1500]);
    
    
    lineatzero(gca,0,[0.1 0.1 0.1]);
    %set(ah_r(s),'Xtick',[-500 -250 0 250 500 750 1000]);
    set(ah_r(s),'XLim',[-500 1000]);
    axis tight
    % set(ah_r(s),'Xticklabel',[-0.5 -0.25 0 0.25 0.5 0.75 1]);
    set(ah_r(s),'Xticklabel',[]);
    
    applyaxprops(gca);
    
    
    % histograms
    % bin data
    
    d=outdat.Tstamps{s};
    
    
    if isempty (d)
        d=zeros(10,1);
    end
    
    ah_h(s)=axes('position',[histpos{1}(s) histpos{2}(s) histpos{3}(s,:)]);
    
    
    % scale bins to Hz?
    if params.normHz{1}
        temphandles=plot_hist_group(d,[],plotparams.edgv,[0 0 0],0,params.normHz);
    else
        temphandles=plot_hist_group(d,[],plotparams.edgv,[0 0 0],1);
    end
    
    set(ah_h(s),'XTick',[-250 0 250 500 750 1000]);
    axis tight
    
    % add color to bars that are significant
    
    if ERstats.pER(s)<params.crit_p
        
        [mini,ix]=min(abs((plotparams.edgv)-mean(ERstats.pwin(s,:))));
        hold all;
        tempdat=zeros(size(plotparams.edgv,2),1);
        tempdat(ix-1)=temphandles.ydata(ix-1);
        
        bar(temphandles.edgv,tempdat,1,'histc');
        
    end
    
    % save temphandles to get ydata
    histhandles{s}=temphandles;
    axis tight
end % end stimuli
%
% rescale y lim to best response
beststimindex=ERstats.beststim;
ylimall=get(ah_h(beststimindex),'Ylim');
set(ah_h(:),'Ylim',ylimall);
% add waveform plot
if params.plotwaveform
    waveformpos=[0.84 0.4 0.14 0.17];
    axwf=axes('position',waveformpos);
    density_plot(spikes, idx);
    ylabel('microV');
    axh.ax_wf=axwf;
end

% axes and figure properties
applyfigprops(fh);
axh.ax_s=ah_s;
axh.ax_r=ah_r;
axh.ax_h=ah_h;
%
set(axh.ax_h(beststimindex).YLabel, 'String','Hz')
set(axh.ax_h(:),'Xtick',[-250 0 250 500 750 1000]);
set(axh.ax_h(beststimindex).XLabel, 'String','ms');
set(axh.ax_h(:),'Box','off');
% rasters
set(axh.ax_r(:),'YColor','none');
set(axh.ax_r(:),'XColor','none');
set(axh.ax_r(:),'Color','none');
% stimuli
set(axh.ax_s(:),'YColor','none');
set(axh.ax_s(:),'XColor','none');
set(axh.ax_s(:),'Color','none');
%
for j=1:nstimuli
    % add lines
    axes(axh.ax_h(j));
    line([0 0],[0 ylimall(end)],'color',[0.1 0.1 0.1],'Linestyle','--');
    
    % change color of bars
    ca=get(axh.ax_h(j),'Children');
    link=findobj(ca, 'Type','Patch');
    for kk=1:numel(link)
        link(kk).FaceColor=[0.5 0.5 0.5];
        link(kk).EdgeColor=[0.5 0.5 0.5];
    end
end