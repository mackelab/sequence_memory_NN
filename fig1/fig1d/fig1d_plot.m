function  fig1d_plot(datafilepath)
%%
load(datafilepath,'bindat','convdat','sortinfo','expinfo');
% tevents to plot
% indices of stimulus_onsets:
% tevents to plot
% stimulus onsets
m=round(median(expinfo.times.stimon));
xtimes=repmat(m,2,1);
ytimes=repmat([1 size(sortinfo,1)]',1, size(xtimes,2));
p=round(median(expinfo.times.panelon)); % round as there is 60ms jitter as to when exactly panel onset occurs


figure;
lineprops={'Color',[0.2 0.2 0.2],'LineStyle','-','LineWidth',0.3};
% rasters for 4 stim positions
loc_ha=axes('position',[0.23 0.53 0.5 0.38]);
%ha=plotRasters(bindat,'x',plotparams_spk.twindow);
ha=plotRasters(bindat,'x',[1 length(bindat)]);

set(loc_ha,'YColor','none');
set(loc_ha,'XColor','none');
set(loc_ha,'Color','none');
% plot 4 stimulus onsets
line(xtimes,ytimes, lineprops{:});
% plot panelonset
line([p p],[1 size(sortinfo,1)],lineprops{:});
% plot lines separating the trials for the different positions
xval=repmat([1 length(bindat)]',1,3);
yval=repmat(find(diff(sortinfo(:,3))),1,2)';
line(xval,yval,lineprops{:});
curxlim=get(gca,'Xlim');
set(gca,'Xlim',([curxlim(1)+500   curxlim(2)]));
set(gca,'Ytick',[yval(1,:)-yval(1)/2 yval(end)+yval(1)/2]);
set(gca,'Yticklabels',{'1','2','3','4'});
ylabel('position');
applyaxprops(gca);
% convolved PSTH
loc_ha=axes('position',[0.23 0.33 0.5 0.38]);
set(loc_ha,'Color','none');
stimpos=unique(sortinfo(:,3));
col=cell2mat(p_colors('yellow_red')');
tvec=((1:length(bindat)) - 1000)-1;
for i=1:numel(stimpos)
    locdat=convdat(:,sortinfo(:,3)==stimpos(i));
    data=mean(locdat,2);
    lpos(i)=plot(tvec,data,'color',col(i,:),'linewidth',1);
    upper=std(locdat,[],2)./sqrt(size(locdat,2))+data;
    lower=data-std(locdat,[],2)./sqrt(size(locdat,2));
    hold on;
    plot(tvec,upper','color',col(i,:),'linewidth',0.2,'Linestyle','-');
    hold on;
    plot(tvec,lower','color',col(i,:),'linewidth',0.2,'Linestyle','-');
    
end

% add stim/panel_onset lines
vline(m-m(1),'k-');
vline(p-p(1),'k-');
xlabel('ms from first stim onset');
ylabel('Hz');
curxlim=get(gca,'Xlim');
curylim=get(gca,'Ylim');
set(gca,'Xlim',([curxlim(1)+500   tvec(end)])); % only show 500ms before stim onset (not full 1000ms)
set(gca,'Xtick',[0 400 800 1200 4000]);
set(gca,'Xticklabels',{'0','400','800','1200','4000'});
set(gca,'Fontsize',10);
set(loc_ha,'Color','none');
set(loc_ha, 'Box','off');
set(loc_ha,'YColor','k');
set(loc_ha,'Ylim',([-1 100]));
set(loc_ha,'Ytick',[0:20:50]);
applyfigprops(gcf);